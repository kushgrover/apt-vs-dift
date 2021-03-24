package heuristics.search;

import explicit.Distribution;
import explicit.ModelExplorer;
import explicit.SMG;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.search.colours.Black;
import heuristics.search.colours.Grey;
import heuristics.search.colours.ModelManager;
import heuristics.search.colours.OccurrenceCounter;
import heuristics.update.StateUpdate;
import parser.State;
import parser.ast.*;
import prism.PrismException;
import org.json.JSONObject;
import prism.PrismLangException;
import prism.PrismSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;



/**
 * Created by maxi on 09.11.2018.
 * Based on BCC+14 atva paper and KKKW18 CAV paper
 * Abstract class to do reachability in black-box, grey-box and white-box setting (settings are then implemented by children-classes.)
 * We have HeuristicSGWhite, which is KKKW18
 * HeuristicSGGrey, which has qualitative knowledge about transitions (not the numbers, but the exact successors)
 * HeuristicSGBlack, which can only sample.
 * Black contains methods to apply the old DQL from BCC+14, and an improved version of DQL.
 * Improvements are
 *  different detection of BECs/BSCCs
 *  handling of BECs (deflating, not collapsing)
 *  learning probs and calculating L/U, not learning L/U directly
 *
 *  Note: This hack does not include collapsing, since it only makes things more complicated and I have no really fast implementation.
 */
public abstract class HeuristicSG extends Heuristic{
    public HeuristicSG(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException {
        super(mc, su, nextState, pme, min);
    }

    //Parameters
    protected boolean collapseMECs = false;
    protected boolean breakImmediately = false;

    /**
     * first three bits are for colour (see HeuristicFactory).
     * fourth bit (%32) is flag whether additional collapsing shall be used
     * fifth bit (%64) decides whether to breakImmediately or sample for longer
     * @param opts
     */
    public void setOpts(int opts){//every bit tells us about an optimization.
        this.collapseMECs = (opts%32)>15;
        this.breakImmediately = (opts%64)>31;
    }

    //for grey and black
    protected OccurrenceCounter oc;

    //path
    protected final Deque<State> visited = new ArrayDeque<>(); //Path
    //protected final Deque<Integer> actions = new ArrayDeque<>(); //actions used; might be needed for black
    //improvement: do not update precomp if not necessary
    protected Set<State> seen = new HashSet<State>();
    protected boolean seenChanged = false;

    //stuff for black explore; these are set anew at the beginning of each explore for black; not used for grey and white.
    protected Map<State,Integer> breakTries = new HashMap<>();
    protected Map<State,Integer> breakThreshold = new HashMap<>();

    public void explore(State currentState) throws PrismException {
        if(this instanceof Black){
            breakTries.clear();
            breakThreshold.clear();
        }
        visited.add(currentState); //add the first state

        while (!stateUpdate.isTarget(currentState) && !stateUpdate.isZero(currentState)) {
            if (exploreBreak(currentState)) {
                if (seenChanged) { //if we have explored new states and are in a setting where we can see the graph
                    // Construct partialMDP from all the states seen so far and (for white/grey) compute 0s/1s
                    updatePrecomp();
                    seenChanged = false;
                }
                handleBEC();//detect SECs in partial game, deflate them
                break; //exit explore loop
            }
            trialSteps++;
            int bestAction = sampleAction(currentState);
            if(bestAction==-1){
                //If we run into a deadlock, we want to set its value
                if(stateUpdate.getQValue(currentState)==null){
                    pme.queryState(currentState);
                    if(stateUpdate.isTarget(currentState)){
                        stateUpdate.setTarget(currentState,true);
                    }
                    else{
                        stateUpdate.setZero(currentState,true);
                    }
                }
                //only happens if there is a deadlock in pme; in that case, we know we can stop exploring;
                //update will handle checking, whether the deadlock state is a 0 or a 1.
                break;
            }
            State succState = nextState.sample(currentState, bestAction);
            seenChanged |= updateCounters(currentState, bestAction, succState);//for white/grey, this never is true; for black, it is true if new trans was found.
            currentState = succState;
            //actions.add(bestAction);
            visited.addLast(currentState); //visited contains the path we take
            seenChanged |= seen.add(currentState);

//            //debug
//            if(debugAll4()){
//                report_smg();
//                System.exit(42);
//            }
        }
    }

    /**
     * Samples the best action for the current state. Depends on whether it is max or min state.
     * For white, uses stateUpdate; for grey/black, uses what they know so far.
     * @param currentState
     * @return
     * @throws PrismException
     */
    public abstract int sampleAction(State currentState) throws PrismException;
    public abstract void updatePrecomp() throws PrismException;//Computing based on pme for grey and white; for black, based on smapling; prob1/0 only for grey/white; for black, we use it to identify EC-candidates
    public abstract boolean exploreBreak(State s) throws PrismException;
    public abstract void handleBEC() throws PrismException;//KKKW18 for grey and white, sampling long time and then deflate; for black: find SEC according to DHKP16-extension, then deflate
    public abstract void update() throws PrismException;//white: KKKW18, BRTDP. grey, black: BCC+14 improvement
    public abstract void doPartialBVI() throws PrismException;
    public abstract boolean updateCounters(State s, int a, State t) throws PrismException; //updates counters in grey/black and returns true, if a new transition was found; needed for black, to ensure we have correct precomp



    //Model
    protected ModelManager mm = new ModelManager(pme, initialState,min);
    public void setCoalition(Set<Integer> coalition, int numPlayers) {
        mm.setCoalition(coalition,numPlayers);
    }


    //Heuristic methods
    @Override
    protected void heuristicStart() throws PrismException {
        modelCheckingTime = System.currentTimeMillis();
        String colour = (this instanceof Black) ? "Black" : ((this instanceof Grey) ? "Grey" : "White");
        mc.getLog().println("HeuristicSG: Version try0\n" +
                colour +
                "\n======================================\n");
        if (verbose) {
            mc.getLog().println("Coalition: " + mm.smg.getCoalitionSet());
            mc.getLog().flush();
        }
        seen.add(initialState);
        if (this instanceof Black){
            mm.buildPartialModelBlack(seen,oc);
        }
        else {
            mm.buildPartialModel(seen);
        }
    }

    @Override
    protected void heuristicStop() throws PrismException {
        if (verbose) {
            long duration = System.currentTimeMillis() - modelCheckingTime;
            mc.getLog().println();
            mc.getLog().println("HeuristicSG model checking time: " + ((double) duration / 1000) + " secs.");
            printProgress(duration);
            //debug
            //report_smg();
        }

        // For APT vs DIFT
        System.out.print("Exporting best and worst case games ... ");
        try{
            exportBestWorst();
        }
        catch(IOException e){
            e.printStackTrace();
        }
        System.out.println("done");
    }

    @Override
    public boolean isDone() throws PrismException {
        StateValue sv = stateUpdate.getQValue(initialState);
        if (sv != null) {
            double lowerBoundInitialState = sv.getLowerBound();
            double upperBoundInitialState = sv.getUpperBound();
            return upperBoundInitialState - lowerBoundInitialState < stateUpdate.getEpsilon() && doneBvi;
        }
        return false;
    }


    boolean doBVI = false; //for grey and black, the guaranteed thing on the statistical estimates;

    //used for time based reporting; might uncomment that stuff again
    //int prints = 0;
    //int checkHowOften =  10000;
    boolean doneBvi = false;
    //Actual algorithm
    @Override
    public void heuristicStep(State state) throws PrismException {
        trials++;
        //Checks every few trials, whether we have been running for some time; if so, we may want to doBVI, and we can printProgress.
        doneBvi = !(this instanceof Grey);//for grey and black, we need to do BVI; white may finish on its own
//        if (trials % checkHowOften == 0) {
//            long duration = System.currentTimeMillis() - modelCheckingTime;
//
//            if(duration>10000*prints) {
//                prints++;
//                printProgress(duration);
//            }
//        }
//        if(trials == 6000000){
//            report_smg();
//            System.exit(32);
//        }

        explore(state);
        update();
        visited.clear();
        
        doBVI = (trials) % 100000 == 0;
        // TODO Pranav: think about when to trigger BVI
        if(doBVI){
            updatePrecomp();
            //report_smg();
            doPartialBVI();
            doneBvi = true;
            //report_smg();
            //System.exit(43);
//            mc.getLog().println("Let's explore a bit");
//            boolean verb = getVerbose();
//            verbose = false;
//            for (int i = 0; i< checkHowOften; i++){
//                explore(state);
//                visited.clear();
//            }
//            verbose = verb;
//            getAvgConf();
        }
    }




    //Output stuff
    protected long modelCheckingTime = 0;
    public void printProgress(long duration) throws PrismException{
        JSONObject json = new JSONObject();
        json.put("Trials", trials);
        json.put("GlobalTimerSecs", ((double) duration / 1000));
        JSONObject value = new JSONObject();
        value.put("Lower", getInitialStateValue().getLowerBound());
        value.put("Upper", getInitialStateValue().getUpperBound());

        json.put("Value", value);
        json.put("Precision", getInitialStateValue().getUpperBound() - getInitialStateValue().getLowerBound());
        json.put("PartialTransDelta",oc.getTransDelta());

        JSONObject actions = new JSONObject();
        int init = mm.smg.getFirstInitialState();
        for (int i = 0; i<mm.smg.getNumChoices(init); i++){
            actions.put("Action" + i, "[" + getLowerBoundActionValue(init,i) + ";" + getUpperBoundActionValue(init,i) + "]");
        }

        json.put("ActionsOfs0", actions);
        json.put("AvgConf", getAvgConf());

        List<State> statesList = mm.smg.getStatesList();
        int numUnset,num00,num11,num01,numWorking,numClose;
        double avgDist;
        numUnset=0;num00=0;num11=0;num01=0;numWorking=0;numClose=0;avgDist=0;
        for (int i = 0; i< mm.smg.getNumStates(); i++){
            StateValue stateValue = stateUpdate.getQValue(statesList.get(i));
            if (stateValue != null) {
                double diff = stateValue.getUpperBound()-stateValue.getLowerBound();
                avgDist+=diff;
                if(diff<0.001){
                    numClose++;
                }
                if(diff==0 &&stateValue.getUpperBound()==0){num00++;}
                else if (diff==0 &&stateValue.getLowerBound()==1){num11++;                }
                else if(diff==1){num01++;}
                else{numWorking++;}
            }
            else{
                numUnset++;
            }
        }
        avgDist = avgDist/mm.smg.getNumStates();

        JSONObject stateInfo = new JSONObject();
        stateInfo.put("numStates",mm.smg.getNumStates());
        stateInfo.put("numUnset",numUnset);
        stateInfo.put("num00",num00);
        stateInfo.put("num11",num11);
        stateInfo.put("num01",num01);
        stateInfo.put("numWorking",numWorking);
        stateInfo.put("numClose",numClose);
        stateInfo.put("avgDist",avgDist);
        json.put("StateInfos", stateInfo);
        mc.getLog().println("JSON: " + json);
    }

    public double getAvgConf() throws PrismException{
        if(!(this instanceof Grey)){
            return -1;
        }
        double confSum = 0;
        int confCount = 0;
        for(State s : oc.allSeenStates()){
            if(stateUpdate.isZero(s)){
                continue;//not interested in the actions of states that are done; they just increase the avg conf value, because we do not improve them.
            }
            pme.queryState(s);
            for(int a = 0; a < pme.getNumChoices(); a++){
                boolean allSuccZero = true;
                for (int t = 0; t < pme.getNumTransitions(a); t++) {
                    allSuccZero &= stateUpdate.isZero(pme.computeTransitionTarget(a, t));
                }
                if(allSuccZero){
                    continue;//see above; we will not use actions that lead to clear zeros, so their confidence will never drop
                }
                confSum += oc.getConfidenceForSApair(s,a);
                confCount++;
            }
        }
        return (confSum/(1.0*confCount));
    }

    //Debug printing
    public void report_smg() throws PrismException{
        mc.getLog().println("FINAL REPORT");
        printProgress(System.currentTimeMillis() - modelCheckingTime);
        int init = mm.smg.getFirstInitialState();
        BitSet done = new BitSet();
        report_dfs(init,done);
        mc.getLog().flush();
    }

    private void report_dfs(int s, BitSet done) throws PrismException{
        SMG smg = mm.smg;
        StateValue value = stateUpdate.getQValue(smg.getStatesList().get(s));
        if(value!=null) {
            mc.getLog().println("State s" + s + " [" + value.getLowerBound() + ";" + value.getUpperBound() + "]" + " called " + smg.getStatesList().get(s) + " player " + smg.getPlayer(s));
        }
        else{
            mc.getLog().println("State s" + s + " [unset]" + " called " + smg.getStatesList().get(s) + " player " + smg.getPlayer(s));
        }
        for (int a = 0; a<smg.getNumChoices(s); a++){
            State sSt = smg.getStatesList().get(s);
            mc.getLog().println("\tAction " + a + ": [" + getLowerBoundActionValue(s,a) + ";" + getUpperBoundActionValue(s,a) + "] {" + oc.getPairCount(sSt,a) + "}");
            Distribution dist_a_i = smg.getChoice(s,a);
            for (int toState : dist_a_i.keySet()) {
                State tSt = smg.getStatesList().get(toState);
                value = stateUpdate.getQValue(tSt);
                if(value!=null) {
                    mc.getLog().println("\t\t" + dist_a_i.get(toState) + " {" + oc.getLowerProbEstimate(sSt,a,tSt) + ", " + oc.getTripleCount(sSt,a,tSt) + "}\t\ts" + toState + "\t\t[" + value.getLowerBound() + ";" + value.getUpperBound() + "]");
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

    protected abstract double getUpperBoundActionValue(int s, int i) throws PrismException;

    protected abstract double getLowerBoundActionValue(int s, int i) throws PrismException;


    ModulesFile mf = mc.getModulesFile();
    private void exportBestWorst() throws IOException, PrismException {
        String dir = mc.getSettings().getString(PrismSettings.OUTPUT_DIR);
        BufferedWriter bestWriter = new BufferedWriter(new FileWriter(dir + "/best.prism"));
        BufferedWriter worstWriter = new BufferedWriter(new FileWriter(dir + "/worst.prism"));

        writeInit(bestWriter);
        writeInit(worstWriter);

        writeModules(bestWriter, worstWriter);

        writeEnd(bestWriter);
        writeEnd(worstWriter);

        bestWriter.close();
        worstWriter.close();
    }

    private void writeInit(BufferedWriter writer) throws IOException, PrismLangException {
        writer.write("smg\n\n");

        // Constants
        ConstantList constList = mf.getConstantList();
        for(int i=0; i<constList.size(); i++){
            writer.write("const int " + constList.getConstantName(i) + " = " + constList.getConstant(i).evaluateInt() + ";\n");
        }
        writer.newLine();

        // Global vars
        for(int i=0; i<mf.getNumGlobals(); i++){
            Declaration var = mf.getGlobal(i);
            Expression low = ((DeclarationInt) var.getDeclType()).getLow();
            Expression high = ((DeclarationInt) var.getDeclType()).getHigh();
            writer.write("global " + var.getName() + ": [" + low + " .. " + high + "];\n");
        }
        writer.newLine();

        // Players
        for(int i=0; i<mf.getNumPlayers(); i++){
            Player p = mf.getPlayer(i);
            writer.write("player " + p.getName() + "\n\t");
            if(i == 0){
                Iterator<String> it = p.getActions().iterator();
                while(it.hasNext()){
                    String act = it.next();
                    writer.write("[" + act + "]");
                    if(it.hasNext())
                        writer.write(", ");
                }
            }
            else if(i == 1){
                writer.write(p.getModules().get(0));
            }
            writer.write("\nendplayer\n\n");
            p.updateSyncs(mf.getSynchs());
        }
    }

    private void writeModules(BufferedWriter bestWriter, BufferedWriter worstWriter) throws IOException, PrismException {
        bestWriter.write("module " + mf.getModule(0).getName() + "\n");
        worstWriter.write("module " + mf.getModule(0).getName() + "\n");

        List<State> allStates = mm.smg.getStatesList();
        Iterator<State> it = allStates.iterator();
        while(it.hasNext()){
            State s = it.next();
            pme.queryState(s);
            for(int a=0; a<pme.getNumChoices(); a++){
                bestWriter.write("\t");
                worstWriter.write("\t");
                writeTransitions(bestWriter, s, a, true);
                writeTransitions(worstWriter, s, a, false);
            }
        }

        bestWriter.write("endmodule\n\n");
        worstWriter.write("endmodule\n\n");
    }

    private void writeTransitions(BufferedWriter writer, State s, int a, boolean bestCase) throws PrismException, IOException {
        int pairCount = oc.getPairCount(s,a);
        List<State> succStates = new ArrayList<State>();
        List<Double> values = new ArrayList<Double>();
        String action = "";
        if((int) s.varValues[1] == 0){
            Player p = mf.getPlayer(0);
            action = p.getSynchs().get(a);
        }
        double sum = 0;
        List<State> allStates = mm.smg.getStatesList();
        Iterator<State> it = allStates.iterator();
        while(it.hasNext()){
            State t = it.next();
            int tripleCount = oc.getTripleCount(s, a, t);
            if(isTrapTransition(s, action ,t)){
                double value = 0;
                if(pairCount != 0)
                    value = (double) tripleCount / (double) pairCount;
                double confWidth = oc.getConfidenceForSApair(s,a);
                if((int) t.varValues[1] == 1){
                    if(bestCase)
                        value -= confWidth;
                    else
                        value += confWidth;
                }
                else{
                    if(bestCase)
                        value += confWidth;
                    else
                        value -= confWidth;
                }
                value = clamp(value);
                sum += value;
                values.add(value);
                succStates.add(t);
            }
            else if (anyOtherTransition(a,t)){
                sum += 1.0;
                values.add(1.0);
                succStates.add(t);
            }
        }

        writer.write("\t[" + action + "] " + getStateString(s) + " -> ");
        for(int i=0; i<succStates.size(); i++){
            if(i < succStates.size()-1)
                writer.write(values.get(i) + " : " + getStateStringPrime(succStates.get(i)) + " + ");
            else
                writer.write(values.get(i) + " : " + getStateStringPrime(succStates.get(i)) + ";\n");
        }
    }

    private boolean anyOtherTransition(int a, State t) throws PrismException {
        State r = pme.getSucc(a);
        if( (int) r.varValues[0] == (int) t.varValues[0] && (int) r.varValues[1] == (int) t.varValues[1]){
            return true;
        }
        return false;
    }

    private double clamp(double x) {
        if(x > 1.0)
            return 1.0;
        else if(x < 0)
            return 0;
        else
            return x;
    }

    private String getStateString(State s) {
        String state = "";
        for(int i=0; i<mf.getNumVars(); i++){
            if(i < mf.getNumVars()-1)
                state += "(" + mf.getVarName(i) + "=" + s.varValues[i].toString() + ") & ";
            else
                state += "(" + mf.getVarName(i) + "=" + s.varValues[i].toString() + ")";
        }
        return state;
    }

    private String getStateStringPrime(State t) {
        String state = "";
        for(int i=0; i<mf.getNumVars(); i++){
            if(i < mf.getNumVars()-1)
                state += "(" + mf.getVarName(i) + "\'" + "=" + t.varValues[i].toString() + ") & ";
            else
                state += "(" + mf.getVarName(i) + "\'" + "=" + t.varValues[i].toString() + ")";
        }
        return state;
    }

    private boolean isTrapTransition(State s, String action, State t) throws PrismLangException {
        if (! action.equals("trap"))
            return false;
        if ((int) s.varValues[0] == (int) t.varValues[0]) {
            if ((int) s.varValues[1] == 0 && (int) t.varValues[1] == 1) {
                return true;
            }
        }
        if ((int) s.varValues[1] == 0 && (int) t.varValues[1] == 0) {
            int n = mf.getConstantList().getConstant(0).evaluateInt();
            if ((int) t.varValues[0] == n+1) {
                return true;
            }
        }
        return false;
    }

    private void writeEnd(BufferedWriter writer) throws IOException {
        for(int i=0; i<mf.getNumRewardStructs(); i++){
            RewardStruct rew = mf.getRewardStruct(i);
            writer.write("rewards \"" + rew.getName() + "\"\n");
            for(int j=0; j<rew.getNumItems(); j++){
                RewardStructItem item = rew.getRewardStructItem(j);
                writer.write("\t" + item.toString() + "\n");
            }
            writer.write("endrewards\n\n");
        }
    }


}
