package heuristics.search.colours;

import explicit.ModelExplorer;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.update.StateUpdate;
import parser.State;
import prism.Pair;
import prism.PrismException;

import java.util.*;

public class Black extends Grey {

    //Note: No collapsing, since that is horrible to implement.

    //Constructor stuff
    public Black(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min, int statesOA, int availableOA, double delta, double pmin, int postMax) throws PrismException {
        super(mc, su, nextState, pme, min, statesOA, availableOA, delta, pmin, postMax);
    }
    @Override
    public void specialStuff(){ //constructor difference between grey and black
        this.oc = new OccurrenceCounter(availableOA, transDelta);
        setStepsUntilSure();
    }

    //variables for caching the MECs
    boolean upToDateMECs = false;
    List<BitSet> currentMECs = null;

    //EXPLORE is fine in grey


    //UPDATE
    @Override
    public boolean updateCounters(State s, int a, State t) {
        return oc.incrementTriple(s, a, t);
    }

    /**
     * Conservative guess of lower/upper bound that remaining prob could go to.
     * Since we do not know whether there is some successor we have not seen yet, always be most conservative.
     * TODO: Small optimization: If we have param postmax (maximum number of succs) we might use the most conservative succ as soon as we have seen postmax many succs.
     * TODO: Or we could override getEstimatedActionBound to not do all the work of remembering the values for the conservative guess; both things are a possible improvement, but they exclude each other
     */
    @Override
    protected double conservativeGuess(boolean lower, double[] allSuccBounds, int numSuccs){
        return lower ? 0 : 1;
    }



    //EXPLOREBREAK and finding ECs for HANDLEBEC
    private int stepsUntilSure;
    protected void setStepsUntilSure(){
        stepsUntilSure = (int) Math.ceil((Math.log(oc.getTransDelta()) / Math.log(1 - pmin))); //from reordering (1-pmin)^seen < transDelta
    }
    /**
     * This method returns false if we should continue sampling,
     *   and true if we are delta-sure that we are stuck in a SEC and should deflate.
     *
     * We break when seeing a state twice as in white. (Running too long is not used)
     * But then we want to make sure, that we are in an EC; so we compute the EC that the currentState is part of,
     *  (according to current knowledge)
     *  and calculate how long we have to sample until we are delta-sure that it is an EC.
     *  Then, for all states in the EC, we set breakThreshold to that value.
     * Then, we return false so sampling is continued;
     *  but the next time we see the state, we do not compute more, until we have seen it breakThreshold many times. (use breakTries to count)
     *  When that is the case, we check again whether it is an EC sure enough, which now has a good chance to be true.
     *  If it is so, we return true; if it is not, we set breakThreshold and breakTries again and continue sampling.
     * If we realize s cannot be in any EC, we set breakTries to -1 and never break for this state anymore.
     */
    @Override
    public boolean exploreBreak(State s) throws PrismException{
        int neededSteps;
        BitSet ec;

        //If we have considered this state before
        if(breakThreshold.containsKey(s)){
            int bt = breakThreshold.get(s);
            if(bt ==-1){
                //can't be in any EC with current strat
                return false;
            }
            //increment counter
            breakTries.replace(s,breakTries.get(s)+1);

            //if counter larger then earlier calculated threshold, check again if we are sure
            if(breakTries.get(s) > breakThreshold.get(s)) {
                Pair<BitSet,Integer> ecNsteps= stepsUntilInECsureEnough(s);
                ec = ecNsteps.first;
                neededSteps = ecNsteps.second;
            }
            //else continue sampling
            else{
                return false;
            }
        }
        //if we have not considered this state before
        else{
            //if we see a cycle
            if(visited.contains(s)){
                //calculate how many steps until we are sure
                Pair<BitSet,Integer> ecNsteps= stepsUntilInECsureEnough(s);
                ec = ecNsteps.first;
                neededSteps = ecNsteps.second;
            }
            //otherwise we have have no reason to do anything, since the explored path is not cycling
            else{
                return false; //continue sampling
            }
        }

        //if we come here, neededSteps and ec have been initialized; all other cases have returned
        //if neededSteps==0, we know that s is in an EC sure enough
        if(neededSteps == 0){
            return true; //is EC sure enough
        }
        //otherwise, we need to update the breakThreshold and breakTries map to allow counting and checking
        //We do this for every state in ec,
        else {
            //if s is not in any SEC, set breakThreshold to -1 to always immediately return false when we see it
            if(ec==null){
                myPut(breakThreshold,s,100);
                myPut(breakTries,s,0); //breakThreshold says how many more visits we need, so we need to restart counting.
            }
            //otherwise, set breakThreshold for all states in ec to neededSteps; since we need the least of them to be seen very often, and the others shouldn't force us to recalculate the ec
            //also initialize breakTries if necessary
            else {
                List<State> int2state = mm.smg.getStatesList();
                for (int i = ec.nextSetBit(0); i >= 0; i = ec.nextSetBit(i + 1)) {
                    myPut(breakThreshold, int2state.get(i), neededSteps); //set threshold
                    myPut(breakTries,int2state.get(i),0); //breakThreshold says how many more visits we need, so we need to restart counting.
                }
            }
            return false; //need to sample more
        }
    }

    /**
     * handles differentiating between replace and put. I don't care, but java does.
     */
    private void myPut(Map<State,Integer> map, State key, int value){
        if(map.containsKey(key)){
            map.replace(key, value);
        }
        else {
            map.put(key, value);
        }
    }

    /**
     * Calculates the SEC containing s according to current strat (this will always be the same for one trial).
     * If it is in no SEC, returns null,-1.
     * If it looks like it is in a SEC, but we are not delta-sure for all actions, return the SEC and the maximum neededSteps, i.e. the difference between stepsUntilSure and pairCount
     *   This is then an underapproximation of how often we have to see the state with that action, because the state might pick a different action.
     *   However, we need to see it at least that often until we need to compute the stepsUntilInECsureEnough again.
     * If it is a SEC sure enough, return the SEC and 0.
     */
    private Pair<BitSet,Integer> stepsUntilInECsureEnough(State s) throws PrismException{
        if (seenChanged) { //if we have explored new states or transitions, need to update the graph to allow sensible SEC computation
            // Construct partialMDP from all the states seen so far
            updatePrecomp();
            seenChanged = false;
        }

        //compute SEC containing s
        BitSet relevantSEC = getSECcontainingState(s);
        //if s is not in any SEC, return (null,100) to indicate that s can never be in a SEC for the current strategy
        // AND with the current view of the model; this is why we do not use sth like -1 to say "never look at this again",
        // but we use 100 to say: If you visit this 100 more times, you are allowed to check whether there is new information
        // and now it actually looks like an EC
        if (relevantSEC == null){
            return new Pair<>(null,100);
        }

        //else we now have the relevantSEC and need to find out how long we need to sample inside this until we are sure
        //for this we find the maximum difference between stepsUntilSure and pairCount over all states in the SEC and staying actions
        int minPairCount = -1;
        //for each state in the SEC
        for(int i = relevantSEC.nextSetBit(0); i >= 0; i = relevantSEC.nextSetBit(i + 1)) {
            State st = mm.smg.getStatesList().get(i);
            //for all staying actions
            for(int a : mm.getAllStayingEC(relevantSEC,i)){
                //get pair count
                int currPairCount = oc.getPairCount(st,a);
                if(minPairCount == -1 || currPairCount < minPairCount){
                    minPairCount = currPairCount; //remember minimal one
                }
            }
        }

        //now we can calculate the maximum difference over all staying state-action pairs
        //We cap it at 0, so that 0 indicates we have seen everything often enough;
        //the only negative value possible for neededSteps is -1 if there is no SEC (see the other return above)
        int neededSteps = Math.max(stepsUntilSure-minPairCount,0);
        if(neededSteps > 0 && neededSteps<10){
            neededSteps=10; //otherwise for small values there is a high probability of not seeing the relevant state-action pair, because we just cycle somewhere else
        }
        return new Pair<>(relevantSEC,neededSteps);
    }

    /**
     * Returns the SEC that contains s, or null if s is not in any SEC according to current strat.
     * Works as usual by computing MECs, and then finding SECs in MEC using the lower bounds.
     * Additionally restricts to that MEC that contains s.
     */
    private BitSet getSECcontainingState(State s) throws PrismException{
        BitSet result = null;

        //compute MECs; sadly, we cannot just compute the MEC containing s
        List<BitSet> mecs = getCurrentMECs();

        //need the lower bounds if we have to calculate SECs
        double[] lowerBounds=null;
        if(!mm.isMDP){
            lowerBounds = mm.boundVectorFromStateUpdate(true,stateUpdate);
        }


        //get number of State s that we need to be in there
        int stateNum = mm.state2Index.get(s);

        for (int i = 0; i < mecs.size(); i++) {
            BitSet mec = mecs.get(i);
            if (!mec.get(stateNum)) {
                continue; //if state not in mec, look at next mec
            }

            //We have found the MEC containing s

            //Now we find the SEC in it that contains s
            if (mm.isMDP) { //in MDP, every MEC is a SEC
                result = mec;
            } else {
                explicit.ECComputer ec = explicit.ECComputer.createECComputer(mc, mm.smg);
                List<BitSet> simpleBCECs = ec.getSimBCECStates(mec, lowerBounds);
                for (int j = 0; j < simpleBCECs.size(); j++) {
                    BitSet sec = simpleBCECs.get(j);
                    if (sec.get(stateNum)) {
                        result = sec;
                    }
                }
            }
        }
        return result;
    }

    /**
     * Handles caching the MECs, because otherwise we waste all our time computing MECs.
     * @return the MECs in the current smg
     */
    private List<BitSet> getCurrentMECs() throws PrismException{
        if(upToDateMECs){
            return currentMECs;
        }

        explicit.ECComputer ec = explicit.ECComputer.createECComputer(mc, mm.smg);

        ec.computeMECStates();
        currentMECs = ec.getMECStates();
        upToDateMECs = true;
        return currentMECs;
    }

    @Override
    public void handleBEC() throws PrismException{
        Pair<BitSet,Integer> ecNsteps= stepsUntilInECsureEnough(visited.getLast());
        if(ecNsteps.second!=0){
            throw new PrismException("HandleBEC called, but we are not delta sure about the SEC. Probably an error in exploreBreak in Black.");
        }
        BitSet SEC = ecNsteps.first;
        //we are delta-sure that we have found a SEC that needs to be deflated.
        //We can use the deflate from white, since that only uses mm.smg, which is something that black sets, too, and that only contains information it may contain.
        //The white deflate also calls stateUpdate.setZero if necessary, so that explore knows this from now on.
        deflate(SEC);
    }


    //UPDATE PRECOMP, or rather partial model; needed for calculating SECs.
    /**
     * Needed for checking whether we have delta-sure EC and for deflate.
     * Should construct the smg, but only include successors which we have seen.
     * This means, that for unplayed actions there are no successors, or rather a dummy successor that shows we have not been there.
     */
    @Override
    public void updatePrecomp() throws PrismException{
        mm.buildPartialModelBlack(seen,oc);
        upToDateMECs = false;
    }

    /**
     * TODO: Better structure, since this is duplicate code from grey
     */
    @Override
    public void doPartialBVI() throws PrismException{
        oc.updateAllProbs(oc.allSeenStates());
        //mc.getLog().println("Running precomputation...");
        updatePrecomp();

        //we can use a better transDelta for partialBVI, since we have less states and know the number of actions for some states
        //in black, differently from grey, we have to assume there are postMax successors for each transition
        int numTrans = 0;
        //get transDelta for current partial model
        for (int s = 0; s<mm.smg.getStatesList().size(); s++){
            for (int a = 0; a<mm.smg.getNumChoices(s); a++){
                numTrans += postMax;
            }
        }
        oc.setTransDelta(delta/numTrans);
        //now we need less steps until we are sure about an EC; we tell this to PartialBVI
        setStepsUntilSure();
        PartialBVI.setStepsUntilSure(stepsUntilSure);
        PartialBVI.progressEpsilon=stateUpdate.getEpsilon();

        //Note: unvisited contains the rim states after we buildPartialModel
        PartialBVI.partialBVI(oc,mm.unvisited,stateUpdate,mm.smg,min,mc.getLog(), true);
        printProgress(System.currentTimeMillis() - modelCheckingTime);

        oc.setTransDelta(transDelta); //set back to old value, whatever that was before
        setStepsUntilSure(); //also set this back to the old value
    }
}
