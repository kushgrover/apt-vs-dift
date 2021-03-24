package heuristics.search.colours;

import explicit.*;
import heuristics.search.StateValue;
import heuristics.update.StateUpdate;
import parser.State;
import prism.PrismException;
import prism.PrismLog;

import java.util.*;

import static java.lang.Math.max;

public class PartialBVI {

    private static SMG smg;
    public  static double progressEpsilon;
    private static OccurrenceCounter oc;
    private static boolean beBlack;
    private static int stepsUntilSure = -1;
    /**
     * Note: Collapsing does not work here, since for games we can't be sure about complex BECs. Also it is difficult to implement a performant version even for MDPs.
     * Hence only deflate, but rim states are not deflated.
     * @param givenoc Used to get probabilities
     * @param rim These states have not been explored, hence they have to stay at [0,1]
     * @param stateUpdate To use "isTarget" and "isZero" and set the initial values
     * @param givenSMG The smg, used for finding BECs and deflating; saved in a static parameter, since all other methods use it
     */
    static void partialBVI(OccurrenceCounter givenoc, BitSet rim, StateUpdate stateUpdate, SMG givenSMG, boolean min, PrismLog log, boolean isBlack) throws PrismException {
        oc = givenoc;
        smg = givenSMG;
        beBlack = isBlack;
        double lowerBounds[], lowerBounds2[], upperBounds[], upperBounds2[];
        int n = smg.getNumStates();
        lowerBounds = new double[n];
        lowerBounds2 = new double[n];
        upperBounds = new double[n];
        upperBounds2 = new double[n];

        //get information from stateUpdate and make BitSets accordingly; also find out unknown states where we need to do stuff
        // Determine set of states actually need to compute values for
        // also initialize upper and lower bounds
        BitSet unknown;
        unknown = new BitSet();
        for (int i = 0; i < n; i++) {
            if (rim.get(i)) { //For black, the only rim state is dummy; this is ensured by ModelManager.buildPartialModelBlack
                //We may not look at this state neither in update nor in deflate; hence not in unknown
                //We only use its bounds to show our absolute ignorance of its value
                lowerBounds[i] = lowerBounds2[i] = 0;
                upperBounds[i] = upperBounds2[i] = 1;
                continue;
            }
            if (stateUpdate.isTarget(index2State(i))) {
                lowerBounds[i] = lowerBounds2[i] = upperBounds[i] = upperBounds2[i] = 1;
            } else if (stateUpdate.isZero(index2State(i))) {
                lowerBounds[i] = lowerBounds2[i] = upperBounds[i] = upperBounds2[i] = 0;
            } else {
                unknown.set(i);
                lowerBounds[i] = lowerBounds2[i] = 0;
                upperBounds[i] = upperBounds2[i] = 1;
            }
        }
        int initialState = smg.getFirstInitialState();
        boolean done = false;
        int iters = 0;

        //Do BVI, using probs from oc, occasionally deflating all things but rim states
        while (!done) {
            iters++;

            bellmanUpdateWithEstimates(lowerBounds, unknown, min, true);
            bellmanUpdateWithEstimates(upperBounds, unknown, min, false);

            //log.println("Average (lowerBounds2, upperBounds2) = \t\t\t(" +
            //        (average(lowerBounds2)) + ", " + (average(upperBounds2)) + ")");
            //log.println("Updated average (lowerBounds, upperBounds) = \t(" +
            //        (average(lowerBounds)) + ", " + (average(upperBounds)) + ")");

            // Check termination
            done = upperBounds[initialState] - lowerBounds[initialState] < progressEpsilon;

            if (iters % 10 == 0) { //arbitrary improvement, not adjust every step, cause it usually takes very long and needs to be propagated
                BitSet copyOfUnknown = new BitSet();
                copyOfUnknown.or(unknown);
                deflate(upperBounds, lowerBounds, copyOfUnknown, min);

                //Also check whether nothing much is happening
                done = !doWeMakeProgress(lowerBounds, lowerBounds2, upperBounds, upperBounds2, initialState);
                //log.println("Before deflate (lowerBounds2, upperBounds2) = \t(" +
                //        (average(lowerBounds2)) + ", " + (average(upperBounds2)) + ")");
                //log.println("After deflate (lowerBounds, upperBounds) = \t\t(" +
                //        (average(lowerBounds)) + ", " + (average(upperBounds)) + ")");
                //Update these two to current estimate, so we can check for progress again next time
                lowerBounds2 = Arrays.copyOf(lowerBounds, lowerBounds.length);
                upperBounds2 = Arrays.copyOf(upperBounds, lowerBounds.length);

                //log.println(iters + " Steps: " + "[" + lowerBounds[initialState] + ";" + upperBounds[initialState] + "]" + " Avg. Lower Bound: " + average(lowerBounds) + " Avg. Upper Bound: " + average(upperBounds));
                //log.flush();
            }
        }

        //Write results into StateValues
        writeResults(stateUpdate, lowerBounds, upperBounds);
    }

    private static State index2State(int i){
        return smg.getStatesList().get(i);
    }

    /**
     * for all states that are not target or zero
     *      bestBound = min ? 1 : 0
     *      for all actions of that state
     *          currentBound = 0
     *          currentProb = 0
     *          worstBound = min ? 1 : 0
     *          for all successor of that action //Careful: Succ might not have been seen
     *              currentBound += bounds[succ] * getConfidenceEstimate(delta, lower)
     *              currentProb += getConfidenceEstimate(delta, lower)
     *              worstBound = (bounds[succ] is worse) ? bounds[succ] : worstBound
     *          currentBound += (1-currentProb) * worstBound; //remaining probability that we are not confident about is taken into account as worst case (but better than just 0/1)
     *          bestBound = (currentBound is better) ? currentBound : bestBound
     *      bounds[state] = bestBound;
     */
    private static void bellmanUpdateWithEstimates(double[] bounds, BitSet unknown, boolean min, boolean lower) throws PrismException{
        for (int s = unknown.nextSetBit(0); s != -1; s = unknown.nextSetBit(s + 1)) {
                double bestBound = getBestBound(s,bounds,min,lower);
                bounds[s] = bestBound;
            }
    }

    private static double getBestBound(int s, double[] bounds, boolean min, boolean lower){
        boolean inCoal = (smg.getPlayer(s)==1);
        boolean stateIsMinimizing = (min ? inCoal : !inCoal);
        double bestBound = stateIsMinimizing ? 1 : 0; //initialize bestBound as worst possible value, so it is updated definitely
        for (int a = 0; a < smg.getNumChoices(s); a++) {
            double currentBound = getBestBound(s,a,bounds,lower);
            bestBound = stateIsMinimizing ? Math.min(bestBound, currentBound) : Math.max(bestBound, currentBound);
        }
        return bestBound;
    }

    private static double getBestBound(int s, int a, double[] bounds, boolean lower){
        double conservativeBound = lower ? 1 : 0; //initialized as the least conservative bounds; then going over the actions takes min/max and returns most conservative bound
        double currentBound = 0;
        double currentProb = 0;
        Iterator<Integer> succsIt = smg.getSuccessorsIterator(s, a);
        while (succsIt.hasNext()) {
            int t = succsIt.next();
            // Using lower prob estimate is safe
            double prob = oc.getLowerProbEstimate(index2State(s), a, index2State(t));

            //TODO: There was the idea of also using avgProb to find out how far we have converged. Use this yet?
            //double prob = oc.getAvgProb(index2State(s), a, index2State(t));
            // if (prob > 0) System.out.println("Greater than 0 prob from oc.getLowerProbEstiamte " + prob);
            // if (prob == 0) System.out.println("Equal to 0 prob from oc.getLowerProbEstiamte " + prob);
            assert 0.0 <= prob && prob <= 1.0;
            double bound = bounds[t];

            currentBound += prob * bound;
            currentProb += prob;
            conservativeBound = lower ? Math.min(conservativeBound, bound) : Math.max(conservativeBound, bound);
        }
        // TODO same possibilities for optimization as in Black.conservativeGuess
        if(beBlack){
            conservativeBound = lower ? 0 : 1;
        }
        currentBound += (1 - currentProb) * conservativeBound;
        return currentBound;
    }

    private static void writeResults(StateUpdate su, double[] l, double[] u){
        for (int i = 0; i < l.length; i++) {
            su.setQValue(index2State(i), new StateValue(l[i], u[i]));
        }
    }

    private static void deflate(double[] u, double[] l, BitSet unknown, boolean min) throws PrismException {
        //compute MECs
        explicit.ECComputerDefault ec = (ECComputerDefault) ECComputer.createECComputer(null, smg);
        ec.computeMECStates(unknown);
        List<BitSet> mecs = ec.getMECStates();

        //for each MEC
        for (BitSet mec : mecs) {
            //Find simBCECs in the MEC
            List<BitSet> simpleBCECs = ec.getSimBCECStates(mec, l);
            for (BitSet simBCEC : simpleBCECs) {
                //In black, we need to ensure that this simBCEC really is an EC, i.e. we have seen every transition often enough
                if(beBlack && !ecSureEnough(simBCEC)){
                    continue; //if we are black and it's not an EC sure enough, we do not deflate. Else we continue as in grey.
                }


                //Note: Need not do anything about rim states, because those are not in unknown
                double bestLeavingUpperBound = getBestLeavingValue(simBCEC, u, min);

                //set all upper bounds to the best upper bound
                for (int s = simBCEC.nextSetBit(0); s >= 0; s = simBCEC.nextSetBit(s + 1)) {
                    double formerValue = u[s];
                    if (formerValue > bestLeavingUpperBound) { //avoids the math min, cause we only enter body, if formerValue is too large
                        u[s] = bestLeavingUpperBound;
                        //System.out.println("BAM!");
                    }
                }
            }
        }
    }

    private static double getBestLeavingValue(BitSet ec,double[] vector, boolean min){
        double bestUpperBoundSoFar = 0;
        //find best outgoing upper bound belonging to player 1
        for (int s = ec.nextSetBit(0); s >= 0; s = ec.nextSetBit(s+1)) {
            boolean inCoal = (smg.getPlayer(s)==1);
            if (min ? !inCoal : inCoal) {
                //mainLog.println("Searching for best leaving value; state belongs to maximizer");
                for (int i = 0; i < smg.getNumChoices(s); i++) {
                    boolean all = smg.allSuccessorsInSet(s, i, ec);
                    //mainLog.println("Action " + i + " all succ in set? " + all);
                    if (!all) {
                        double upperBound = 0;
                        for (int succ : smg.getChoice(s,i).keySet()){
                            upperBound +=  getBestBound(s,i,vector,false)* vector[succ]; //only called for deflating, hence always upper bound, hence lower is false
                        }
                        if (upperBound>bestUpperBoundSoFar){
                            bestUpperBoundSoFar = upperBound;
                        }
                    }
                }
            }
        }
        //Target should never be given to this method, since it is not in the unknown BitSet
//        if (bestUpperBoundSoFar==0){
//            //Check for target in simBCEC
//            for (int s = ec.nextSetBit(0); s >= 0; s = ec.nextSetBit(s+1)) {
//                if (vector[s]==1) {
//                    bestUpperBoundSoFar=1;//if we find target in simBCEC, all states in there should have value 1.
//                    break;
//                }
//            }
//        }
        return bestUpperBoundSoFar;
    }

    private static boolean ecSureEnough(BitSet candidate){
        int minPairCount = -1;
        //for each state in the SEC
        for(int i = candidate.nextSetBit(0); i >= 0; i = candidate.nextSetBit(i + 1)) {
            State st = smg.getStatesList().get(i);
            //for all staying actions
            for(int a : getAllStayingEC(candidate,i)){
                //get pair count
                int currPairCount = oc.getPairCount(st,a);
                if(minPairCount == -1 || currPairCount < minPairCount){
                    minPairCount = currPairCount; //remember minimal one
                }
            }
        }

        return minPairCount>=stepsUntilSure;
    }

    /**
     * Copy pasted from ModelManager. TODO: Make the whole PartialBVI structure nicer; probably have a static util class with all methods that are used twice.
     */
    public static Set<Integer> getAllStayingEC(BitSet relevantSEC, int j) {
        Set<Integer> result = new HashSet<>();
        for(int a = 0; a < smg.getNumChoices(j);a++) {
            boolean all = smg.allSuccessorsInSet(j, a, relevantSEC);
            if(all){
                result.add(a);
            }
        }
        return result;
    }

    public static void setStepsUntilSure(int steps){
        stepsUntilSure=steps;
    }

    //Looks whether there is some progress in this BVI; needed, because if rim states are important or confidence is low, we need a different break criterion than convergence
    //bounds iters to 1/progressEpsilon * how often this is called * numStates or sth like that
    private static boolean doWeMakeProgress(double[] lowerBounds, double[] lowerBounds2, double[] upperBounds, double[] upperBounds2, int initialState){
        if (lowerBounds[initialState] - lowerBounds2[initialState] > progressEpsilon || upperBounds2[initialState] - upperBounds[initialState] > progressEpsilon) {
            return true; //easy, even in startstate we see progress;
        }
        if (average(lowerBounds) - average(lowerBounds2) > progressEpsilon || average(upperBounds2) - average(upperBounds) > progressEpsilon) {
            return true; //one of the bound functions is making progress in average (i.e. somewhere, just needs to be propagated); values have to be higher here, since this is aggregated
        }
        return false;
    }

    public static double sum(double[] vec) {
        double result = 0;
        for (double v : vec) {
            result += v;
        }
        return result;
    }

    public static double average(double[] vec) {
        return sum(vec) / vec.length;
    }

}
