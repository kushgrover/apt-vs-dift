package heuristics.search;

import explicit.*;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.update.StateUpdate;
import parser.State;
import parser.ast.Expression;
import prism.Prism;
import prism.PrismException;
import heuristics.CachedModelExplorer;
//import heuristics.CollapseModelExplorer;

import java.util.*;

/**
 * Created by maxi on 5/28/17.
 * Based on the code by mateusz-ujma for his thesis "On Verification and Controller Synthesis for Probabilistic Systems at Runtime"
 */
public class HeuristicRTDP_Adj extends Heuristic{

    private long modelCheckingTime = 0;
    private Set<State> seen = new HashSet<State>();
    private SMG smg = new SMG();
    public void setCoalition(Set<Integer> coalition,int numPlayers){
        if(coalition.isEmpty()){ //for mer
            coalition.add(1);
        }
        if(numPlayers==0){//for mdps
            numPlayers=1;
            coalition.add(1);
        }
        smg.setCoalitionSet(coalition);
        smg.setCoalitionPlayerMapFromIntSet(coalition,numPlayers);
    }
    private Map<State, Integer> state2Index = new HashMap<State, Integer>();
    private BitSet unvisited = new BitSet();
    private HashSet<State> zeros = new HashSet<>();

    private boolean seenChanged = false;
    private boolean adjustedSth = true;
    private BitSet finished = new BitSet();
    private final Deque<State> visited = new ArrayDeque<>();
    private int exploreSteps = 0;
    private int adjustCount = 0;
    private int adjustThreshold = 0;
    private static final int PROGRESS_REPORT_TRIALS = 100;
    public void setOpts(int opts){//every bit tells us about an optimization. TODO: Explain all here. (also in prism.java)
        this.repeatedAdjustment = (opts%2)>0;//if the least bit is one, we do repeated adjustment
        this.collapseMECs = (opts%4)>1;
        this.breakImmediately = (opts%8)>3;
        this.onlyLastEc = (opts%16)>7;
        this.k_i_mat = (opts%32)>15;
        this.k_i_large = (opts%64)>31;
        adjustThreshold = ((opts%128)>63) ? 1000 : 0;
        mc.getLog().println("repeatedAdjustment: " + repeatedAdjustment + " collapseMECS: " + collapseMECs + " breakImmediately: " + breakImmediately + " onlyLastEc: " + onlyLastEc
                + " k_i_mat: " + k_i_mat + " k_i_large: " + k_i_large + " adjustmentThreshold: " + ((opts%128)>63));
        mc.getLog().flush();
    }
    private boolean repeatedAdjustment = false; //set all opts to false at beginning; they are calculated once the opts-int is set by setOpts
    private boolean collapseMECs = false;
    private boolean breakImmediately = false;
    private boolean onlyLastEc = false;
    private boolean k_i_mat=false;
    private boolean k_i_large=false;
    private boolean isMDP=false;
    public void setIsMDP(boolean isMDP){
        this.isMDP=isMDP;
    }
//    private boolean breakEarly=false;
//    private static boolean noticedProblem = false;
//    public static void noticeProblem(){noticedProblem=true;}//when we preferLeaving and realize that a non-leaving action is preferred to all leaving actions, we know we have to adjust
//    private HashMap<State,BitSet> state2leavingActions = new HashMap<State,BitSet>();

    public HeuristicRTDP_Adj(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException {
        super(mc, su, nextState, pme, min);
    }

    @Override
    protected void heuristicStart() throws PrismException{
        modelCheckingTime = System.currentTimeMillis();
        mc.getLog().println("RTDP_Adj Version: Second Bench Try for CAV");
        if (verbose) {
            mc.getLog().println("Coalition: " + smg.getCoalitionSet());
            mc.getLog().flush();
        }
        smg = buildPartialModel();
    }

    @Override
    protected void heuristicStop() throws PrismException{
        long duration = System.currentTimeMillis() - modelCheckingTime;
        if (verbose) {
            mc.getLog().println();
            mc.getLog().println("Heuristic RTDP_Adj model checking time in " + ((double) duration / 1000) + " secs.");
            //debug
            //report_smg();
        }
    }


    //old comment; TODO: look into this when I'm done with basics <= important now! what happens with min, max and minmax? Also: Bounded Properties; look at prism_examples
    //min is never used; but most probably it should be; look what happens if pctl contains min, max or minmax
    //looks like p1 is always maxing and p2 always minning, or sth; maxmin and minmax dont change behaviour
    //however, since updatePrecomp does use min, we can achieve nontermination, since updatePrecomp says were done (have 0 node), but update(s) picks a maximizing action
    //21.06: taken from tobi
    //21.06.: Tobis code (and Mateusz') break explore the first time they see a state twice (to update more often?). They only do collapse every 1000 steps
    //  should we rather explore, until we see one state 5 (or 1000) times, and then break explore, adjust, back to explore?
    //30.06. Pseudocode suggest we run til we see target, zero, or we have a certain amount of trial steps
    @Override
    public void heuristicStep(State state) throws PrismException {
        trials++;

//        if(trials>500000){
//            report_smg();
//            System.exit(3);
//        }

//        if(trials==100){
//            boolean debugBreak = true;
//            report_smg();
//        }

        if (trials%100==0 && verbose){
            long duration = System.currentTimeMillis() - modelCheckingTime;
            mc.getLog().println("Trials: " + trials + " Global timer: " + ((double) duration / 1000) + " secs.");
            printProgress(smg,stateUpdate);
            //mc.getLog().println(debugOutput);
            //System.exit(3);
        }

        // EXPLORE Phase
        explore(state);

        // UPDATE Phase
        update();
    }

    //private String debugOutput;

    public void explore(State currentState) throws PrismException{
        while (!stateUpdate.isTarget(currentState) && !stateUpdate.isZero(currentState) && !exploreBreak(currentState)) {
            seenChanged |= seen.add(currentState);
//            mc.getLog().print((state2Index.containsKey(currentState)?state2Index.get(currentState):"n")+",");
//            mc.getLog().flush();
            visited.addLast(currentState); //visited contains the path we take
            int bestAction=-1;
//            if (preferLeaving) {
//                bestAction = stateUpdate.update(currentState, state2leavingActions);
//            } else {
                //TODO: If I ever look at this again, need to adjust this to stateUpdate changes bestAction = stateUpdate.update(currentState);
                //debugOutput += Integer.toString(bestAction);
//            }
            currentState = nextState.sample(currentState, bestAction);

            trialSteps++;
            exploreSteps++;
        }
//        mc.getLog().println();

        //debugOutput += "\n";
        //debugPrint(currentState);

        if (exploreBreak(currentState)) { //if we ran in cycle (so there is some EC that possibly needs adjusting) or every 100 steps, so we do it even if we reach target all the time

//            if(preferLeaving){
//                noticedProblem=false;
//            }
            if (seenChanged) { //if we have explored new states
                // Construct partialMDP from all the states seen so far
                updatePrecomp();
                seenChanged = false;
            }
            //debugPrint(currentState);

            //remember how many trials since last adjust
            adjustCount++;
            //adjust probabilities if number exceeds threshold; new heuristic (15.01.18): Let Threshold depend on L
            if(adjustCount>adjustThreshold) {
                if(onlyLastEc){
                    adjustProbabilities(smg,findLastEc(visited));
                }
                adjustProbabilities(smg);
                adjustCount = 0;
            }
            //empty visited, so we can jump to explore again
            visited.removeAll(new ArrayDeque<>(visited));
            exploreSteps = 0;

            //mc.getLog().print("Current SMG: ");
            //printQValues(smg,stateUpdate);
            if (verbose) {
                printProgress(smg, stateUpdate);
            }

        }
        else { //we have reached a terminal
            stateUpdate.setTargetOrZeroValue(currentState);
//            if (verbose) {
//                if (stateUpdate.isTarget(currentState)) {
//                    mc.getLog().println("Reached target state " + currentState);
//                }
//                if (stateUpdate.isZero(currentState)) {
//                    mc.getLog().println("Reached zero state " + currentState);
//                    zeros.add(currentState);
//                }
//            }
        }
    }

    public boolean exploreBreak(State currentState) throws PrismException{
        if(breakImmediately){
            return visited.contains(currentState) && visited.size()>30 ;// || (preferLeaving&&noticedProblem); //TODO: delete hack (the min path length); was only for cloud;
        }
//        if(breakEarly){
//            StateValue qval = stateUpdate.getQValue(currentState);
//            if (qval!= null && qval.getUpperBound()-qval.getLowerBound() ==0){
//                return true;
//            }
//        }
        if(k_i_mat){
            if(k_i_large){return exploreSteps> (10^(5 + ((int)(trials/10))));}
            else{return exploreSteps> (10^(1 + ((int)(trials/10))))*smg.getNumStates();}
        }
        else{
            if(k_i_large){return exploreSteps> (500*smg.getNumStates());}// || (preferLeaving&&noticedProblem);}
            else{return exploreSteps> (2*smg.getNumStates());}// || (preferLeaving&&noticedProblem);}
        }

    }

    private BitSet findLastEc(Deque<State> visited){
        //Checks for the largest EC that we have in visited; if there is none, we adjust nothing, else we find the largest simBCEC that we are stuck in
        BitSet result=new BitSet();

        BitSet seen = new BitSet();

        int lastState = state2Index.get(visited.removeLast());

        while(!visited.isEmpty()) {
            //iterate over visited from behind
            int s = state2Index.get(visited.removeLast());
            seen.set(s);
            if (result.get(s) || s == lastState) { //first run: Find EC once we find the last state again; after that: Find anything from result
                result.or(seen);
            }
        }
        return result;

    }

    public void debugPrint(State currentState) throws PrismException{
        if (exploreBreak(currentState)){
            mc.getLog().println("Broken");
        }
        else{
            mc.getLog().println("Terminal");
        }
        for (State s : visited){
            StateValue value = stateUpdate.getQValue(s);
            if(value!=null) {
                mc.getLog().println("State s" + s + " [" + value.getLowerBound() + ";" + value.getUpperBound() + "]" + " player " + stateUpdate.getPlayer(s)); //player acc to stateUpdate
            }
            else{
                mc.getLog().println("State s" + s + " [unset]"  + " player (acc to SU)" + stateUpdate.getPlayer(s));
            }
        }
        mc.getLog().println("State s" + currentState + " Target: " + stateUpdate.isTarget(currentState));
    }


    public void update() throws PrismException{
        while (!visited.isEmpty()) {
            State visitedState = visited.removeLast();
            // Upper and lower bounds are propagated from the successors
            //TODO: If I ever look at this again, need to adjust this to stateUpdate changes stateUpdate.update(visitedState);
        }

        if ((trials + 1) % PROGRESS_REPORT_TRIALS == 0) {
            reportProgress(trials, trialSteps);
        }
        exploreSteps = 0;

//        if(trials>10000){
//            report_smg();
//            System.exit(3);
//        }
//        if (isDone()) {
//            report_smg();
//            SMGModelChecker smgModelChecker = getMC();
//            smg = buildPartialModel();
//            //smg.findDeadlocks(true);
//            BitSet target = computeTarget(smgModelChecker, smg);
//            HashSet<State> tar = new HashSet<State>();
//            for (int i=target.nextSetBit(0);i>=0;i = target.nextSetBit(i+1)){
//                tar.add(smg.getStatesList().get(i));
//            }
//            mc.getLog().println(tar);
//            mc.getLog().println("\nResult: " + sv.getUpperBound());
//            heuristicStop();
//            System.exit(3);
//        }
    }

    public void report_smg() throws PrismException{
        mc.getLog().println("FINAL REPORT");
        printProgress(smg,stateUpdate);
        int init = smg.getFirstInitialState();
        BitSet done = new BitSet();
        report_dfs(init,done);
        mc.getLog().flush();
    }

    private void report_dfs(int s, BitSet done) throws PrismException{//TODO: add zero and target check
        StateValue value = stateUpdate.getQValue(smg.getStatesList().get(s));
        if(value!=null) {
            mc.getLog().println("State s" + s + " [" + value.getLowerBound() + ";" + value.getUpperBound() + "]" + " called " + smg.getStatesList().get(s) + " player " + smg.getPlayer(s));
        }
        else{
            mc.getLog().println("State s" + s + " [unset]" + " called " + smg.getStatesList().get(s) + " player " + smg.getPlayer(s));
        }
        for (int i = 0; i<smg.getNumChoices(s); i++){
            mc.getLog().println("\tAction " + i + ": [" + stateUpdate.getLowerBoundActionValue(smg.getStatesList().get(s),i).getLowerBound() + ";" + stateUpdate.getUpperBoundActionValue(smg.getStatesList().get(s),i).getUpperBound() + "]");
            Distribution dist_a_i = smg.getChoice(s,i);
            for (int toState : dist_a_i.keySet()) {
                value = stateUpdate.getQValue(smg.getStatesList().get(toState));
                if(value!=null) {
                    mc.getLog().println("\t\t" + dist_a_i.get(toState) + "\t\ts" + toState + "\t\t[" + value.getLowerBound() + ";" + value.getUpperBound() + "]");
                }
                else{
                    mc.getLog().println("\t\t" + dist_a_i.get(toState) + "\t\ts" + toState + "\t\t[unset]");
                }
            }
        }
        mc.getLog().println();
        done.set(s);

        for (int i = 0; i<smg.getNumChoices(s); i++){
            Distribution dist_a_i = smg.getChoice(s,i);
            for (int toState : dist_a_i.keySet()) {
                if(!done.get(toState)) {
                    report_dfs(toState, done);
                }
            }
        }
    }

    @Override
    public boolean isDone() throws PrismException {
        StateValue sv = stateUpdate.getQValue(initialState);
        if(sv != null) {
            double lowerBoundInitialState = sv.getLowerBound();
            double upperBoundInitialState = sv.getUpperBound();
            return upperBoundInitialState - lowerBoundInitialState < stateUpdate.getEpsilon();
        }
        return false;
    }

    private SMGModelChecker smgModelChecker = getMC();
    //from pranav
    private void updatePrecomp() throws PrismException {
        // Search for MECs in the partial MDP and collapse them
        if (collapseMECs) {
            collapseMECs(smg, unvisited);
        }
        smg = buildPartialModel();
        smg.findDeadlocks(true);
        BitSet target = computeTarget(smgModelChecker, smg);
        //if (!target.isEmpty()) {
            BitSet prob0Target = new BitSet();
            prob0Target.or(target);
            prob0Target.or(unvisited);

            BitSet prob0 = computeProb0(smgModelChecker, smg, prob0Target);
            for (int i = prob0.nextSetBit(0); i >= 0; i = prob0.nextSetBit(i + 1)) {
                State s = smg.getStatesList().get(i);
                stateUpdate.setZero(s, true);
                finished.set(i);
            }
            BitSet prob1 = computeProb1(smgModelChecker, smg, target);
            for (int i = prob1.nextSetBit(0); i >= 0; i = prob1.nextSetBit(i + 1)) {
                State s = smg.getStatesList().get(i);
                if (seen.contains(s)) {
                    stateUpdate.setTarget(s, true);
                    finished.set(i);
                }
            }
        //}
    }



    private SMG buildPartialModel() throws PrismException{
        List<State> statesList = smg.getStatesList();
        if(statesList == null) {
            statesList = new ArrayList<State>();
        }
        unvisited.clear();
        if(smg.getFirstInitialState() == -1) {
            //Empty model
            int index = smg.addState();
            smg.addInitialState(index);
            state2Index.put(initialState, index);
            statesList.add(initialState);

            pme.queryState(initialState);
            smg.setPlayer(index, pme.getPlayerForState());
        }
        Iterator<State> it = seen.iterator();
        while(it.hasNext()) {
            State s = it.next();
            Integer index = state2Index.get(s);
            if(index == null) {
                index = smg.addState();
                state2Index.put(s, index);
                statesList.add(s);

                pme.queryState(s);
                smg.setPlayer(index, pme.getPlayerForState());
            }
            List<Distribution> dists = buildAllDistributions(smg, s, statesList);
            smg.clearState(index);
            for(int i=0;i<dists.size();i++) {
                smg.addChoice(index, dists.get(i));
            }
        }
        smg.setStatesList(statesList);
        if (verbose) {
            mc.getLog().println("Model built with " + smg.getNumStates() + " states");
            long duration = System.currentTimeMillis() - modelCheckingTime;
            mc.getLog().println("Global timer: " + ((double) duration / 1000) + " secs.");
            mc.getLog().flush();
        }
        return smg;
    }


    //created this, since I want my partial model to have all transitions for all states available; to check, whether there is one leaving the MEC
    private List<Distribution> buildAllDistributions(SMG smg, State s, List<State> statesList) throws PrismException{
        List<Distribution> dists = new ArrayList<Distribution>();
        pme.queryState(s);
        int choices = pme.getNumChoices();
        for(int i=0;i<choices;i++) {
            pme.queryState(s);
            Distribution d = new Distribution();
            int trans = pme.getNumTransitions(i);
            for(int j=0;j<trans;j++) {
                pme.queryState(s);
                double prob = pme.getTransitionProbability(i, j);
                State t = pme.computeTransitionTarget(i,j);
                Integer index = state2Index.get(t);
                if(index == null) {
                    index = smg.addState();
                    state2Index.put(t, index);
                    statesList.add(t);

                    pme.queryState(t);
                    smg.setPlayer(index, pme.getPlayerForState());
                }
                d.add(index, prob);
                if(!seen.contains(t)) {
                    unvisited.set(index);
                }
            }
            dists.add(d);
        }
        return dists;
    }

    private BitSet computeProb0(SMGModelChecker smgModelChecker, SMG smg, BitSet target) throws PrismException{
        return smgModelChecker.prob0(smg, null, target, min, !min);
    }

    private BitSet computeProb1(SMGModelChecker smgModelChecker, SMG smg, BitSet target) throws PrismException{
        return smgModelChecker.prob1(smg, null, target, min, !min);
    }

    private BitSet computeTarget(SMGModelChecker smgModelChecker, SMG smg) throws PrismException {
        Expression targetExp = null;
        if(mc.getExpression().getOperand1() != null) {
            targetExp = mc.getExpression();
        } else {
            targetExp = mc.getExpression().getOperand2();
        }
        BitSet statesOfInterest = new BitSet(smg.getNumStates());
        statesOfInterest.set(0,smg.getNumStates());
        return smgModelChecker.checkExpression(smg, targetExp,statesOfInterest).getBitSet();
    }


    private SMGModelChecker getMC() throws PrismException{
        SMGModelChecker smgModelChecker = new SMGModelChecker(mc);
        smgModelChecker.setModulesFileAndPropertiesFile(mc.getModulesFile(), mc.getPropertiesFile());
        return smgModelChecker;
    }


    private void adjustProbabilities(SMG smg) throws PrismException {
        if (verbose) {
            mc.getLog().println("Starting adjusting probabilities");
            mc.getLog().flush();
        }
        long start = System.currentTimeMillis();

        //get BitSet for all states that are not finished yet; if a state is finished, prob0 or prob1 will take care of propagating that to all, that are in a MEC with it
        BitSet all = new BitSet();
        //I don't see why we need to ignore the ones which are not one (upper bound < 1.0); probably only in the case of 1 player ECs or sth like that.
//        BitSet notOne = new BitSet();
        for(int i=0;i<smg.getNumStates();i++) {
            all.set(i);
//            StateValue sv = stateUpdate.getQValue(smg.getStatesList().get(i));
//            if(sv != null) {
//                if(sv.getUpperBound() < 1.0) {
//                    notOne.set(i);
//                }
//            }
        }
//        all.xor(notOne);

        all.xor(finished);

        //compute MECs
        explicit.ECComputer ec = explicit.ECComputer.createECComputer(mc, smg);

        ec.computeMECStates(all);
        List<BitSet> mecs = ec.getMECStates();

        //need the lower bounds if we have to calculate simBCECs
        double[] lowerBounds=null;
        if(!isMDP){
            //full lowerBounds needed for each simBCEC calculation
            lowerBounds = boundVectorFromStateUpdate(true);
        }

        //loop at least once
        boolean firstTime = true;
        boolean doRepeat  = repeatedAdjustment;//only do more than one repetition if the opts number tells us to.
        while(firstTime || (adjustedSth && doRepeat)) {
            firstTime = false;
            adjustedSth = false;

            for (int i = 0; i < mecs.size(); i++) {
                BitSet mec = mecs.get(i);

                //if deadlock, ignore; if bottom circle controlled CEC, prob0 will handle it; checking for isTarget or isZero or isDeadlock explicitly takes too long; so just mix the cases and ignore them
                if (mec.cardinality() == 1) { //we should only ignore trap states. We may not ignore single loops along the path, those are important. If it is not left, trap state or not explored => COntinue. If there are actions leaving it: adjust
                    List<Map<State, Double>> actions = getAllLeavingMEC(smg, mec, true);
                    if(actions.size()==0){
                        continue;
                    }
                }



                if(isMDP) {//every MEC is a simBCEC in an MDP
                    adjustProbabilities(smg, mec);
                }
                else{
                    List<BitSet> simpleBCECs = ec.getSimBCECStates(mec, lowerBounds);
                    for (int j = 0; j<simpleBCECs.size();j++){
                        BitSet simBCEC = simpleBCECs.get(j);
                        adjustProbabilities(smg, simBCEC);
                    }
                }
            }

        }

        long duration = System.currentTimeMillis() - start;
        if (verbose) {
            mc.getLog().println("Adjusting probabilities done in " + (double) duration / 1000 + " secs.");
            mc.getLog().flush();
        }

    }

    private double[] boundVectorFromStateUpdate(boolean L){
        double[] result = new double[smg.getNumStates()];
        for (int i = 0; i<smg.getNumStates();i++){
            StateValue val = stateUpdate.getQValue(smg.getStatesList().get(i));
            if (val!=null) {
                if(L) {
                    result[i] = stateUpdate.getQValue(smg.getStatesList().get(i)).getLowerBound();
                }
                else{
                    result[i] = stateUpdate.getQValue(smg.getStatesList().get(i)).getUpperBound();
                }
            }
            else{
                if(L) {
                    result[i] = 0; //know nothing yet
                }
                else{
                    result[i] = 1;
                }
            }
        }
        return result;
    }

    private boolean containsTarget(SMG smg, BitSet mec) throws PrismException{
        for (int i = mec.nextSetBit(0); i >= 0; i = mec.nextSetBit(i+1)) {
            State s = smg.getStatesList().get(i);
            if(stateUpdate.isTarget(s)) {
                return true;
            }
        }
        return false;
    }

    private boolean containsZero(SMG smg, BitSet mec) throws PrismException{
        for (int i = mec.nextSetBit(0); i >= 0; i = mec.nextSetBit(i+1)) {
            State s = smg.getStatesList().get(i);
            if(stateUpdate.isZero(s)) {
                return true;
            }
        }
        return false;
    }


    private void adjustProbabilities(SMG smg, BitSet mec) throws PrismException{

        double bestUpperBoundSoFar = 0;
        //find best outgoing upper bound belonging to player 1
        for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)) {
            if (smg.getPlayer(s)==1) {
                for (int i = 0; i < smg.getNumChoices(s); i++) {
                    boolean all = smg.allSuccessorsInSet(s, i, mec);
                    if (!all) {
                        double upperBound = stateUpdate.getUpperBoundActionValue(smg.getStatesList().get(s),i).getUpperBound();
                        if (upperBound>bestUpperBoundSoFar){
                            bestUpperBoundSoFar = upperBound;
                        }
                    }
                }
            }
        }

        //set all upper bounds to the best upper bound
        for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)) {
            StateValue formerValue = stateUpdate.getQValue(smg.getStatesList().get(s));
            double formerU = (formerValue==null) ? 1 : formerValue.getUpperBound();
            if(formerU>bestUpperBoundSoFar) { //avoids the math min, cause we only enter body, if formerValue is too large
                stateUpdate.setQValue(smg.getStatesList().get(s), new StateValue(formerValue.getLowerBound(), Math.min(formerValue.getUpperBound(), bestUpperBoundSoFar)));
                if(formerU - bestUpperBoundSoFar > stateUpdate.getEpsilon()){

                    adjustedSth = true;
                }
                if (verbose) {
                    mc.getLog().println("Bam! Useful adjustment!");// on state " + smg.getStatesList().get(s) + " former value " + formerValue.getUpperBound() + " new upper bound " + bestUpperBoundSoFar);
                }
            }
        }

    }

    //I only care about actions of player 1 leaving (or, more abstractly, actions of the players in the coalition); that is what onlyCoalition is for
    //If it is true, I only return actions by players in the coalition
    private List<Map<State, Double>> getAllLeavingMEC(SMG smg, BitSet mec, boolean onlyCoalition) throws PrismException {
        List<Map<State, Double>> actions = new ArrayList<Map<State, Double>>();
        for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)) {
            for(int i = 0; i < smg.getNumChoices(s);i++) {
                if(onlyCoalition && !(smg.getPlayer(s)==1)){
                    //we only want states from players in the coalition AND this state does not belong to the coalition
                    continue;
                }
                boolean all = smg.allSuccessorsInSet(s, i, mec);
                if(!all) {
                    actions.add(getDistribution(smg, s, i));
//                    if(preferLeaving && !state2leavingActions.containsKey(smg.getStatesList().get(s))){ //add the leaving action for this state, so we rather choose it when selecting actions in explore; this only has to be done one time
//                        addToLeavingActions(smg.getStatesList().get(s),i);
//                    }
                }
            }
        }
        return actions;
    }

//    private void addToLeavingActions(State s, int action){
//        if (state2leavingActions.containsKey(s)){
//            state2leavingActions.get(s).set(action);
//        }
//        else{
//            BitSet newSet = new BitSet();
//            newSet.set(action);
//            state2leavingActions.put(s,newSet);
//        }
//    }

    private Map<State, Double> getDistribution(SMG smg, int s, int i) {
        Map<State, Double> d = new HashMap<State, Double>();
        Iterator<Map.Entry<Integer,Double>> it =  smg.getTransitionsIterator(s, i);
        while(it.hasNext()) {
            Map.Entry<Integer, Double> e = it.next();
            d.put(smg.getStatesList().get(e.getKey()), e.getValue());
        }
        return d;
    }

    private void printQValues(SMG smg, StateUpdate stateupdate){
        List<State> statesList = smg.getStatesList();
        for (int i = 0; i< smg.getNumStates(); i++){
            StateValue stateValue = stateUpdate.getQValue(statesList.get(i));
            if (stateValue != null) {
                mc.getLog().print("State " + statesList.get(i) + ": [" + stateValue.getLowerBound() + "|" + stateValue.getUpperBound() + "]\t");
            }
            else{
                mc.getLog().print("State " + statesList.get(i) + ": Unset\t");
            }
        }
        mc.getLog().println();
    }

    private void printProgress(SMG smg, StateUpdate stateUpdate) throws PrismException{
        List<State> statesList = smg.getStatesList();
        int numUnset,num00,num11,num01,numWorking,numClose;
        double avgDist;
        numUnset=0;num00=0;num11=0;num01=0;numWorking=0;numClose=0;avgDist=0;
        for (int i = 0; i< smg.getNumStates(); i++){
            StateValue stateValue = stateUpdate.getQValue(statesList.get(i));
            if (stateValue != null) {
                double diff = stateValue.getUpperBound()-stateValue.getLowerBound();
                avgDist+=diff;
                if(diff<0.001){numClose++;}
                if(diff==0 &&stateValue.getUpperBound()==0){num00++;}
                else if (diff==0 &&stateValue.getLowerBound()==1){num11++;
                    /*int choices = smg.getNumChoices(i);
                    boolean everythingOk = false;
                    for (int j=0; i<choices;i++){
                        everythingOk = everythingOk || (stateUpdate.getLowerBoundActionValue(statesList.get(i), j).getLowerBound()==stateValue.getLowerBound());
                    }
                    if (!everythingOk){
                        boolean isTarget = stateUpdate.isTarget(statesList.get(i));
                        int k = 1;
                    }*/
                }
                else if(diff==1){num01++;}
                else{numWorking++;}
            }
            else{
                numUnset++;
            }
        }
        avgDist = avgDist/smg.getNumStates();
        mc.getLog().println("numStates: " + smg.getNumStates() + "\tnumUnset: " + numUnset + "\tnum00: " + num00 + "\tnum11: " + num11 + "\tnum01: " + num01 + "\tnumWorking: " + numWorking + "\tnumClose: " + numClose + "\tavgDist: " + avgDist);
        mc.getLog().println("Actions of s0: ");
        int init = smg.getFirstInitialState();
        for (int i = 0; i<smg.getNumChoices(init); i++){
            mc.getLog().println("\tAction " + i + ": [" + stateUpdate.getLowerBoundActionValue(smg.getStatesList().get(init),i).getLowerBound() + ";" + stateUpdate.getUpperBoundActionValue(smg.getStatesList().get(init),i).getUpperBound() + "]");
        }
//        mc.getLog().print("Finished: ");
//        for (int i = finished.nextSetBit(0); i>=0; i = finished.nextSetBit(i+1)){
//            mc.getLog().print(statesList.get(i));
//        }
//        mc.getLog().println();
        mc.getLog().flush();
    }

    private void collapseMECs(SMG smg, BitSet unvisited) throws PrismException {
        if (verbose) {
            mc.getLog().println("Starting Simple BCEC collapsing...");
            mc.getLog().flush();
        }
        long start = System.currentTimeMillis();
        BitSet all = new BitSet();
        BitSet notOne = new BitSet();
        for(int i=0;i<smg.getNumStates();i++) {
            all.set(i);
            StateValue sv = stateUpdate.getQValue(smg.getStatesList().get(i));
            if(sv != null) {
                if(sv.getUpperBound() < 1.0) {
                    notOne.set(i);
                }
            }
        }
        //all.xor(notOne);
        all.xor(unvisited);
        explicit.ECComputerDefault ec = (ECComputerDefault) ECComputer.createECComputer(mc, smg);
        ec.computeMECStates(all);

        //full lowerBounds needed for each simBCEC calculation
        double[] lowerBounds = boundVectorFromStateUpdate(true);
        double[] upperBounds = boundVectorFromStateUpdate(false);

        List<BitSet> mecs = ec.getMECStates();
        for(int i=0;i<mecs.size();i++) {
            BitSet mec = mecs.get(i);
            if(isMDP){
                List<Map<State, Double>> actions = getAllLeavingMEC(smg, mec,false);
                collapse(smg,mec,actions);
            }
            else {
                List<BitSet> simpleBCECs = ec.getSimBCECStates(mec, lowerBounds);
                for (int j = 0; j < simpleBCECs.size(); j++) {
                    BitSet simBCEC = simpleBCECs.get(j);
                    if (surelySimple(smg, simBCEC, upperBounds, lowerBounds, ec)) {
                        List<Map<State, Double>> actions = getAllLeavingMEC(smg, simBCEC, false);
                        collapse(smg, mec, actions);
                    }
                }
            }

        }
        long duration = System.currentTimeMillis() - start;
        if(verbose) {
            mc.getLog().println("Simple BCEC collapsing done " + (double) duration / 1000 + " secs.");
            mc.getLog().flush();
        }
    }

    private boolean surelySimple(SMG smg, BitSet simBCEC,double[] upperBounds,double[] lowerBounds,explicit.ECComputerDefault ecC) throws PrismException{
        boolean result = true;
        //for all states of player circ, all staying actions must have lower U than the highest L of leaving ones
        for (int s = simBCEC.nextSetBit(0); s >= 0 && result; s = simBCEC.nextSetBit(s+1)) {
            if(!(stateUpdate.getPlayer(smg.getStatesList().get(s))==1)){
                //minimizers state, need to check
                double minStayingU = ecC.getMinStayingValue(s,simBCEC,smg,upperBounds);
                double maxLeavingL = ecC.getMaxLeavingValue(s,simBCEC,smg,lowerBounds);
                result = result && minStayingU <= maxLeavingL;
            }
        }

        return result;
    }


    private void collapse(SMG smg, BitSet mec, List<Map<State, Double>> actions) throws PrismException {
        if (pme instanceof CachedModelExplorer) {
            for (int j = mec.nextSetBit(0); j >= 0; j = mec.nextSetBit(j + 1)) {
                collapse(smg, j, actions);
            }
        }
//        else if (pme instanceof CollapseModelExplorer) {
//            ((CollapseModelExplorer) pme).collapse(smg, mec, actions);
//        }
    }

    private void collapse(SMG smg, int s, List<Map<State, Double>> actions) throws PrismException {
        State state = smg.getStatesList().get(s);
        if (pme instanceof CachedModelExplorer) {
            CachedModelExplorer cme = (CachedModelExplorer)pme;
            cme.updateCache(state, actions);
        }
    }


}