package heuristics.search.colours;

import explicit.ModelExplorer;
import parser.State;
import prism.PrismException;

import java.util.*;

/**
 * Used by grey and black to count the number we have seen a state action pair (s,a), and the number we have seen a transition to target state t.
 * //TODO: We might switch everything from states to ints; also in the simulation/update/precomp, if possible. Use ints everywhere, and only transform to states if necessary
 */
public class OccurrenceCounter {

  /**
     * From https://stackoverflow.com/questions/81346/most-efficient-way-to-increment-a-map-value-in-java
     * Noted there as fastest method to increment value in a HashMap.
     */
    class MutableInt {
        int value = 1; // note that we start at 1 since we're counting
        public void increment () { ++value;      }
        public int  get ()       { return value; }
    }

    private int availableOA;
    private double transDelta;
    private ModelExplorer pme=null;

    OccurrenceCounter(int availableOA, double transDelta){
        this.availableOA=availableOA;
        this.transDelta =transDelta;
    }

    OccurrenceCounter(int availableOA, double transDelta, ModelExplorer pme){
        this.availableOA=availableOA;
        this.transDelta =transDelta;
        this.pme=pme;
    }

    // stores the number of times each s, a pair has been encountered
    private Map<State, int[]> pairCounter = new HashMap<>();
    // stores the number of times each s, a, s triple has been encountered
    private Map<State, List<Map<State,MutableInt>>> tripleCounter = new HashMap<>();
    //We also save the probabilities in the occurrenceCounter, since it has everything it needs to compute them right there
    private Map<State, List<Map<State,Double>>> lowerProbEstimates = new HashMap<>();
    private Map<State, double[]> remainingProb = new HashMap<>();

    /**
     * Increments occurrences of triple (s,a,t) as well as the ones of pair (s,a).
     * Does not touch probabilities.
     * Note: If we collapse, stuff needs to adapted; BUT: for black collapsing really is very hard, since in EC we need to play until we reach exit; in grey, it is kinda possible, but still weird...
     * Returns true, if new triple was created (new trans found, seenChanged, need to updatePrecomp)
     */
    boolean incrementTriple(State s, int a, State t){
        if (pairCounter.containsKey(s)) { //then tripleCounter contains s as well
            pairCounter.get(s)[a] = pairCounter.get(s)[a]+1;

            if (tripleCounter.get(s).get(a) == null) {
                tripleCounter.get(s).set(a, makeTargetCounter(t));
                return true; //knew state, but not action
            } else {
                if (tripleCounter.get(s).get(a).get(t) == null) {
                    tripleCounter.get(s).get(a).put(t, new MutableInt());
                    return true; //knew action, but not trans
                } else {
                    tripleCounter.get(s).get(a).get(t).increment();
                    return false; //knew whole triple
                }
            }
        } else { //make new list and map for pair and triple counters
            int[] pairValues = new int[availableOA];
            for (int i = 0; i < pairValues.length; i++) {
                pairValues[i] = 0;
            }
            pairValues[a] = 1;
            pairCounter.put(s, pairValues);

            List<Map<State, MutableInt>> tripleValues = new ArrayList<>(availableOA);
            initializeArrayList(tripleValues);
            tripleValues.set(a, makeTargetCounter(t));
            tripleCounter.put(s, tripleValues);
            return true; //didn't know state
        }
    }

    //helper method for incrementTriple, because all things might not be initialized yet; handles one more indirection
    //returns a map from target state to occurrences. Needs to be inserted as tripleCounter.get(s).get(a).put(t,this)
    private Map<State,MutableInt> makeTargetCounter(State t){
        Map<State,MutableInt> targetCount = new HashMap<>();
        targetCount.put(t,new MutableInt());
        return targetCount;
    }
    //Needed, because otherwise I can't add something at position 1 before adding something at 0.
    private <A,B> void  initializeArrayList(List<Map<A,B>> l){
        for (int i=0; i<availableOA;i++){
            l.add(i,null);
        }
    }

    /**
     * Updates all the probability estimates with the latest counter values.
     * This is sensible, because after changing a single triple, all probs for this state-action pair change (since the total number has changed)
     * Currently, we have enough memory, so we can set and store all lowerProbEstimates to access them faster, since we might access same prob very often. Some sensible caching method would be good in the long run
     * We give this a stateset to update, which is exactly the path we have seen. Then by this we always update all the states where some prob has changed, and everything accessed is initialized (since we only access stuff that is in the path)
     */
    void updateAllProbs(Set<State> path) throws PrismException{
        for (State s : path){//earlier was tripleCounter.keySet() //path


            // Pranav: for each seen action
            for (int a = 0; a < tripleCounter.get(s).size(); a++){
                if (tripleCounter.get(s).get(a) == null) {
                    continue; //action not played; only exists since a larger action has been played before; leave at null to indicate this.
                }
                double confidenceWidth = getConfidenceWidth(pairCounter.get(s)[a]);
                double probSum = 0;
                // Pranav: for each successor which has been seen for an action (s, a)
                for (State t : tripleCounter.get(s).get(a).keySet()) { //t certainly seen before, since in keyset
                    probSum += updateSingleProb(s, a, t, confidenceWidth);
                }
                if(remainingProb.containsKey(s)) {
                    remainingProb.get(s)[a] = 1 - probSum; //After adding the sum for all successors, I know which probability remains.
                }
                else{
                    double[] actionValues = new double[availableOA];
                    Arrays.fill(actionValues, 1.0); //remainingProb for an action I have not played is 1.
                    actionValues[a] = 1 - probSum;
                    remainingProb.put(s,actionValues);
                }
            }
        }
    }


    /**
     * Like incrementTriple, but it only updates the lowerProbEstimates and leaves the counters alone
     * @return lowerProbEstimate for (s, a, t), using confidenceWidth.
     */
    private double updateSingleProb(State s, int a, State t, double confidenceWidth) throws PrismException{
        //Calculate the lower estimate
        double prob = -1;//needs to be initialized, because else the compiler whines. Actually is never used.
        boolean doCalculate = true;
        if(pme!=null){//in grey setting, use qualitative information
            pme.queryState(s);
            if (pme.getNumTransitions(a)==1){
                prob = 1; //this assumes we call this with the correct t. This assumption is safe, since we only call it, if s,a,t has been observed.
                doCalculate = false;
            }
        }
        if(doCalculate){ //only becomes false if we are in grey and have a trivial case
            prob = Math.max(0, ((double) tripleCounter.get(s).get(a).get(t).get())/pairCounter.get(s)[a] - confidenceWidth); //Average - confidenceWidth; capped to not go below 0.
        }

        //put lower estimate into the cache
        if(lowerProbEstimates.containsKey(s)){
            if(lowerProbEstimates.get(s).get(a)!=null){
                lowerProbEstimates.get(s).get(a).put(t,prob);
            }
            else{
                lowerProbEstimates.get(s).set(a, makeProbStore(t, prob));
            }
        }
        else{
            List<Map<State, Double>> tripleValues = new ArrayList<>();
            initializeArrayList(tripleValues);
            tripleValues.set(a, makeProbStore(t,prob));
            lowerProbEstimates.put(s, tripleValues);
        }
        return prob;
    }

    /**
     * Utility, again removing an indirection
     */
    private Map<State,Double> makeProbStore(State t, double prob){
        Map<State,Double> probStore = new HashMap<>();
        probStore.put(t,prob);
        return probStore;
    }
    /**
     * Utility.
     * Returns confidenceWidth according to Hoeffding bound; error prob is transDelta, pairCount is the number of experiments.
     */
    private double getConfidenceWidth(int pairCount){
        return Math.sqrt((Math.log(transDelta))/(-2 * pairCount));
    }

    /**
     * Returns the probability of a single transition.
     * More precisely, our lower estimate of the probability that is correct with 1-transDelta.
     * Handles nulls and everything; only method that other classes can use to get probabilities
     */
    public double getLowerProbEstimate(State s, int a, State t){
        if(lowerProbEstimates.containsKey(s)){
            if(lowerProbEstimates.get(s).get(a)!=null){
                if(lowerProbEstimates.get(s).get(a).containsKey(t)){
                    return lowerProbEstimates.get(s).get(a).get(t);
                }
            }
        }
        return 0; //in all other cases, I have not seen the state or action or target before, so I cannot make any statement about the probability, except that it is at least 0.
    }

    /**
     * Returns the probability that is left over for a state-action pair when adding all lower estimates
     */
    double getRemainingProb(State s, int a) throws PrismException{
        if(remainingProb.containsKey(s)){
            return remainingProb.get(s)[a];
        }
        return 1; //if not initialized, then there are no succs, and hence all prob remains
        //throw new PrismException("Remaining prob accessed without being initialized for state " + s + ". Have you called updateAllProbs?");
    }

    /**
     * Returns all successors we have seen for a state action pair so far. Thus all states that can have a positive probability.
     */
    Set<State> getSeenSuccs(State s, int a){
        if (tripleCounter.containsKey(s)){
            if(tripleCounter.get(s).get(a)!=null){
                return tripleCounter.get(s).get(a).keySet();
            }
        }
        return new HashSet<>();
    }

    public void setTransDelta(double transDelta){
        this.transDelta=transDelta;
    }
    public double getTransDelta(){
        return transDelta;
    }

    public int getPairCount(State s, int a){
        if(pairCounter.containsKey(s)) {
            return pairCounter.get(s)[a];
        }
        else{
            return 0;
        }
    }

    public int getTripleCount(State s, int a, State t){
        if(tripleCounter.keySet().contains(s)){
            if(tripleCounter.get(s).get(a)!=null){
                if(tripleCounter.get(s).get(a).containsKey(t)) {
                    return tripleCounter.get(s).get(a).get(t).value;
                }
            }
        }
        return 0;
    }

    Set<State> getTransitionTargets(State s, int a){
        if(tripleCounter.keySet().contains(s)){
            if(tripleCounter.get(s).get(a)!=null){
                return tripleCounter.get(s).get(a).keySet();
            }
        }
        return new HashSet<State>();
    }
    public Set<State> allSeenStates() {
        return tripleCounter.keySet();
    }

    public double getConfidenceForSApair(State s, int a){
        if(pairCounter.containsKey(s)){
            if(pairCounter.get(s)[a]!=0) {
                return Math.min(getConfidenceWidth(pairCounter.get(s)[a]),1);
            }
        }
        return 1;
    }

    public double getAvgProb(State s, int a, State t){
        if(lowerProbEstimates.containsKey(s)){
            if(lowerProbEstimates.get(s).get(a)!=null){
                if(lowerProbEstimates.get(s).get(a).containsKey(t)){
                    return (1.0 * tripleCounter.get(s).get(a).get(t).value) / pairCounter.get(s)[a];
                }
            }
        }
        return 0;
    }

}
