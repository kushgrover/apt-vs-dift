package heuristics.search.colours;

import explicit.Distribution;
import explicit.ModelExplorer;
import explicit.SMG;
import explicit.SMGModelChecker;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.search.StateValue;
import heuristics.update.StateUpdate;
import parser.State;
import parser.ast.Expression;
import prism.PrismException;

import java.util.*;

public class ModelManager {
    public ModelManager(ModelExplorer pme, State initialState, boolean min) {
        this.pme=pme;
        this.initialState=initialState;
        this.min=min;
    }


    //Handling partial model
    protected BitSet unvisited = new BitSet();
    public Map<State, Integer> state2Index = new HashMap<State, Integer>();
    protected ModelExplorer pme;
    protected State initialState;

    /**
     * Model and coalition stuff
     */
    protected boolean isMDP=false;
    public SMG smg = new SMG(); //safety risk, but makes things easier to read
    boolean min;
    public void setCoalition(Set<Integer> coalition, int numPlayers) {
        if (coalition.isEmpty()) { //for mer
            coalition.add(1);
        }
        if (numPlayers == 0) {//for mdps
            numPlayers = 1;
            coalition.add(1);
            isMDP=true;
        }
        smg.setCoalitionSet(coalition);
        smg.setCoalitionPlayerMapFromIntSet(coalition, numPlayers);
    }

    public boolean isMaxState(int s){
        boolean inCoal = smg.getPlayer(s)==1;
        return min ? !inCoal : inCoal;
    }
    public boolean isMaxState(State s) throws PrismException{
        pme.queryState(s);
        int player = pme.getPlayerForState();
        boolean inCoal = smg.getCoalitionSet().contains(pme.getPlayerForState());
        return min ? !inCoal : inCoal;
    }

    /**
     * Builds the partial model with all information there. It is up to the children-colour classes to ignore the information.
     */
    public SMG buildPartialModel(Set<State> seen) throws PrismException {
        List<State> statesList = smg.getStatesList();
        if (statesList == null) {
            statesList = new ArrayList<State>();
        }
        unvisited.clear();
        if (smg.getFirstInitialState() == -1) {
            //Empty model
            int index = smg.addState();
            smg.addInitialState(index);
            state2Index.put(initialState, index);
            statesList.add(initialState);

            pme.queryState(initialState);
            smg.setPlayer(index, pme.getPlayerForState());
        }
        Iterator<State> it = seen.iterator();
        while (it.hasNext()) {
            State s = it.next();
            Integer index = state2Index.get(s);
            if (index == null) {
                index = smg.addState();
                state2Index.put(s, index);
                statesList.add(s);

                pme.queryState(s);
                smg.setPlayer(index, pme.getPlayerForState());
            }
            List<Distribution> dists = buildAllDistributions(s, statesList, seen);
            smg.clearState(index);
            for (int i = 0; i < dists.size(); i++) {
                smg.addChoice(index, dists.get(i));
            }
        }
        smg.setStatesList(statesList);
        return smg;
    }

    private List<Distribution> buildAllDistributions(State s, List<State> statesList, Set<State> seen) throws PrismException {
        List<Distribution> dists = new ArrayList<Distribution>();
        pme.queryState(s);
        int choices = pme.getNumChoices();
        for (int i = 0; i < choices; i++) {
            pme.queryState(s);
            Distribution d = new Distribution();
            int trans = pme.getNumTransitions(i);
            for (int j = 0; j < trans; j++) {
                pme.queryState(s);
                double prob = pme.getTransitionProbability(i, j);
                State t = pme.computeTransitionTarget(i, j);
                Integer index = state2Index.get(t);
                if (index == null) {
                    index = smg.addState();
                    state2Index.put(t, index);
                    statesList.add(t);

                    pme.queryState(t);
                    smg.setPlayer(index, pme.getPlayerForState());
                }
                d.add(index, prob);
                if (!seen.contains(t)) {
                    unvisited.set(index);
                }
            }
            dists.add(d);
        }
        return dists;
    }

    /**
     * Builds the partial model for the black setting.
     * Restrictions are not to use probabilities (or only estimates, but we do not need them, I think),
     * and to only have successors we have seen before.
     * This implies we have a dummy state as successor of everything I have not seen so far.
     * I argue at relevant positions, why the knowledge we use is ok to use. This looks like this //C: <argument>
     */
    public SMG buildPartialModelBlack(Set<State> seen, OccurrenceCounter oc) throws PrismException {
        if(dummy==-1){
            setDummyState();
            unvisited.clear();
            unvisited.set(dummy);
        }
        assert unvisited.cardinality()==1; //always there should only be dummy in this


        List<State> statesList = smg.getStatesList();
        if (statesList == null) {
            statesList = new ArrayList<State>();
        }
        if (smg.getFirstInitialState() == -1) {
            //Empty model
            int index = smg.addState();
            smg.addInitialState(index);
            state2Index.put(initialState, index);
            statesList.add(initialState);

            pme.queryState(initialState); //C: we know the intial state and its player
            smg.setPlayer(index, pme.getPlayerForState());
        }
        Iterator<State> it = seen.iterator();
        while (it.hasNext()) {
            State s = it.next();
            Integer index = state2Index.get(s);
            if (index == null) {
                index = smg.addState();
                state2Index.put(s, index);
                statesList.add(s);

                pme.queryState(s);//C: For any state we have seen, we know its player
                smg.setPlayer(index, pme.getPlayerForState());
            }
            //Only difference to normal buildPartialModel TODO: If I have time, restructure the methods so this obviously is the only difference
            List<Distribution> dists = buildAllDistributionsBlack(s, statesList, oc);
            smg.clearState(index);
            for (int i = 0; i < dists.size(); i++) {
                smg.addChoice(index, dists.get(i));
            }
        }
        smg.setStatesList(statesList);
        return smg;
    }

    /**
     * For black, there has to be a dummy state where all unseen actions lead. We fix a unique such state in the beginning and then use it the whole remaining time.
     */
    int dummy = -1;
    public void setDummyState(){
        if(dummy!=-1){
            return;
        }
        else{
            //have to fix a dummy state
            dummy = smg.addState();
            smg.getStatesList().add(dummy,new State(0));
            Distribution d = new Distribution();
            d.add(dummy,1);
            smg.addChoice(dummy,d);
        }
    }
    public int getDummyState() throws PrismException{
        if(dummy==-1){
            throw new PrismException("Dummy not set yet!");
        }
        return dummy;
    }


    private List<Distribution> buildAllDistributionsBlack(State s, List<State> statesList, OccurrenceCounter oc) throws PrismException {
        List<Distribution> dists = new ArrayList<Distribution>();
        pme.queryState(s);
        int choices = pme.getNumChoices();//C: We know the available actions
        for (int a = 0; a < choices; a++) {
            Distribution d = new Distribution();

            Set<State> targets = oc.getTransitionTargets(s,a);

            //if we have not observed this action yet, make a transition to dummy
            if(targets.size()==0){
                assert dummy != -1; //dummy has to be set for this.
                d.add(dummy,1);
                continue; //this action is done now, continue to next one
            }

            //otherwise, we add transitions as observed
            for (State t : targets) {
                double prob = 1; //Just saying that this successor certainly exists, saying nothing about probs. LowerProbEstimates or avgProb are also possible
                Integer index = state2Index.get(t);
                if (index == null) {
                    index = smg.addState();
                    state2Index.put(t, index);
                    statesList.add(t);

                    pme.queryState(t);//C: we know the player for all states
                    smg.setPlayer(index, pme.getPlayerForState());
                }
                d.add(index, prob);
            }
            dists.add(d);
        }
        return dists;
    }


    /**
     * Stuff for precomputing targets and sinks
     */
    protected SMGModelChecker smgModelChecker;
    protected BitSet finished = new BitSet();

    public void updatePrecomp(Set<State> seen, HeuristicsSMGModelChecker mc, StateUpdate stateUpdate) throws PrismException {
        setMC(mc);
        buildPartialModel(seen);
        smg.findDeadlocks(true);
        BitSet target = computeTarget(smgModelChecker, smg, mc);
        //if (!target.isEmpty()) {
        BitSet prob0Target = new BitSet();
        prob0Target.or(target);
        prob0Target.or(unvisited);

        BitSet prob0 = computeProb0(smgModelChecker, smg, prob0Target,min);
        for (int i = prob0.nextSetBit(0); i >= 0; i = prob0.nextSetBit(i + 1)) {
            State s = smg.getStatesList().get(i);
            stateUpdate.setZero(s, true);
            finished.set(i);
        }
        BitSet prob1 = computeProb1(smgModelChecker, smg, target,min);
        for (int i = prob1.nextSetBit(0); i >= 0; i = prob1.nextSetBit(i + 1)) {
          State s = smg.getStatesList().get(i);
          if (seen.contains(s)) {
            stateUpdate.setTarget(s, true);
            finished.set(i);
          }
        }
    }
    private void setMC(HeuristicsSMGModelChecker mc) throws PrismException{
        if (smgModelChecker==null) {
            smgModelChecker = new SMGModelChecker(mc);
            smgModelChecker.setModulesFileAndPropertiesFile(mc.getModulesFile(), mc.getPropertiesFile());
        }
    }
    private BitSet computeProb0(SMGModelChecker smgModelChecker, SMG smg, BitSet target, boolean min) throws PrismException{
        return smgModelChecker.prob0(smg, null, target, min, !min);
    }

    private BitSet computeProb1(SMGModelChecker smgModelChecker, SMG smg, BitSet target, boolean min) throws PrismException{
        return smgModelChecker.prob1(smg, null, target, min, !min);
    }

    private BitSet computeTarget(SMGModelChecker smgModelChecker, SMG smg, HeuristicsSMGModelChecker mc) throws PrismException {
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



    /**
     * Stuff for handling BECs
     */

    //I only care about actions of player 1 leaving (or, more abstractly, actions of the players in the coalition); that is what onlyCoalition is for
    //If it is true, I only return actions by players in the coalition
    public List<Map<State, Double>> getAllLeavingMEC(BitSet mec, boolean onlyCoalition) throws PrismException {
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

    private Map<State, Double> getDistribution(SMG smg, int s, int i) {
        Map<State, Double> d = new HashMap<State, Double>();
        Iterator<Map.Entry<Integer,Double>> it =  smg.getTransitionsIterator(s, i);
        while(it.hasNext()) {
            Map.Entry<Integer, Double> e = it.next();
            d.put(smg.getStatesList().get(e.getKey()), e.getValue());
        }
        return d;
    }

    public double[] boundVectorFromStateUpdate(boolean L, StateUpdate stateUpdate){
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

    private boolean surelySimple(SMG smg, BitSet simBCEC,double[] upperBounds,double[] lowerBounds,explicit.ECComputerDefault ecC, StateUpdate stateUpdate, boolean min) throws PrismException{
        boolean result = true;
        //for all states of player circ, all staying actions must have lower U than the highest L of leaving ones
        for (int s = simBCEC.nextSetBit(0); s >= 0 && result; s = simBCEC.nextSetBit(s+1)) {
            boolean inCoalition = stateUpdate.getPlayer(smg.getStatesList().get(s))==1;
            boolean isMinimizerState = min ? inCoalition : !inCoalition;
            if(isMinimizerState){
                double minStayingU = ecC.getMinStayingValue(s,simBCEC,smg,upperBounds);
                double maxLeavingL = ecC.getMaxLeavingValue(s,simBCEC,smg,lowerBounds);
                result = result && minStayingU <= maxLeavingL;
            }
        }

        return result;
    }

    public Set<Integer> getAllStayingEC(BitSet relevantSEC, int j) {
        Set<Integer> result = new HashSet<>();
        for(int a = 0; a < smg.getNumChoices(j);a++) {
            boolean all = smg.allSuccessorsInSet(j, a, relevantSEC);
            if(all){
                result.add(a);
            }
        }
        return result;
    }
}
