package heuristics.search.colours;

import explicit.*;
import heuristics.CachedModelExplorer;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.search.HeuristicSG;
import heuristics.search.StateValue;
import heuristics.update.StateUpdate;
import parser.State;
import parser.ast.Expression;
import prism.PrismException;

import java.util.*;

public class White extends HeuristicSG {

    public White(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException {
        super(mc, su, nextState, pme, min);
        oc = new OccurrenceCounter(-1,-1);
    }

    @Override
    public void update() throws PrismException{
        while (!visited.isEmpty()) {
            State visitedState = visited.removeLast();
            // Upper and lower bounds are propagated from the successors
            stateUpdate.update(visitedState, mm.isMaxState(visitedState));
        }
    }

    @Override
    public boolean exploreBreak(State currentState) throws PrismException{
        if(breakImmediately){
            return visited.contains(currentState);
        }
        else{
            return visited.size()> (2*mm.smg.getNumStates());
        }
    }

    @Override
    public int sampleAction(State currentState) throws PrismException{
        return stateUpdate.sampleBestAction(currentState, mm.isMaxState(currentState));
    }

    @Override
    public void handleBEC() throws PrismException{
        if (verbose) {
            mc.getLog().println("Starting deflating");
            mc.getLog().flush();
        }
        long start = System.currentTimeMillis();

        //get BitSet for all states that are not finished yet; if a state is finished, prob0 or prob1 will take care of propagating that to all, that are in a MEC with it
        BitSet all = new BitSet();
        for(int i=0;i<mm.smg.getNumStates();i++) {
            all.set(i);
        }
        all.xor(mm.finished);

        //compute MECs
        explicit.ECComputer ec = explicit.ECComputer.createECComputer(mc, mm.smg);

        ec.computeMECStates(all);
        List<BitSet> mecs = ec.getMECStates();

        //need the lower bounds if we have to calculate simBCECs
        double[] lowerBounds=null;
        if(!mm.isMDP){
            lowerBounds = mm.boundVectorFromStateUpdate(true,stateUpdate);
        }

        for (int i = 0; i < mecs.size(); i++) {
            BitSet mec = mecs.get(i);

            //if deadlock, ignore; if bottom circle controlled CEC, prob0 will handle it; checking for isTarget or isZero or isDeadlock explicitly takes too long; so just mix the cases and ignore them
            //Note: This is different for black, because there we have no prob0
            if (mec.cardinality() == 1) { //we should only ignore trap states. We may not ignore single loops along the path, those are important. If it is not left, trap state or not explored => COntinue. If there are actions leaving it: adjust
                List<Map<State, Double>> actions = mm.getAllLeavingMEC(mec, true);
                if (actions.size() == 0) {
                    continue;
                }
            }


            if (mm.isMDP) {//every MEC is a simBCEC in an MDP
                deflate(mec);
            } else {
                List<BitSet> simpleBCECs = ec.getSimBCECStates(mec, lowerBounds);
                for (int j = 0; j < simpleBCECs.size(); j++) {
                    BitSet simBCEC = simpleBCECs.get(j);
                    deflate(simBCEC);
                }
            }
        }

        long duration = System.currentTimeMillis() - start;
        if (verbose) {
            mc.getLog().println("SECs done in " + (double) duration / 1000 + " secs.");
            mc.getLog().flush();
        }
    }

    /**
     * Given an EC (it needn't be one, but it is only called with such),
     * this finds the best exit (acc. to U) of maximizer and sets the upper bound of all states in the EC to U(bestExit).
     * An argument why this method is ok for grey and black as well is at all positions where sth unallowed seems to happen, marked with //C:
     */
    protected void deflate(BitSet mec) throws PrismException{
        //find best outgoing upper bound belonging to Maximizer
        double bestUpperBoundSoFar = 0;
        for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)) {
            if (mm.isMaxState(s)) { //C: We always know the player of a state
                for (int i = 0; i < mm.smg.getNumChoices(s); i++) { //C: We always know the number of actions
                    boolean all = mm.smg.allSuccessorsInSet(s, i, mec); //C: In grey, this is fine, since it is qualitative. In black, the smg is constructed to contain only the successors that we have seen so far, so we do not use disallowed knowledge
                    if (!all) {
                        double upperBound = getUpperBoundActionValue(s,i);//C: Own method to calculate this based on estimates in grey/black
                        if (upperBound>bestUpperBoundSoFar){
                            bestUpperBoundSoFar = upperBound;
                        }
                    }
                }
            }
        }
        if(targetInMec(mec)){//C: We know whether a state is target
            bestUpperBoundSoFar=1; //if target, set to 1. TODO: Check this first and if target, don't compute the stuff
        }

        //For black: Our only method of finding sinks is deflate. Hence deflate needs to call stateUpdate.setZero, if the upperBound is deflated to 0.
        boolean setZero = false;
        if(bestUpperBoundSoFar==0){
            setZero = true;
        }
        //set all upper bounds to the best upper bound
        for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)) {
            StateValue formerValue = stateUpdate.getQValue(mm.smg.getStatesList().get(s));
            double formerU = (formerValue==null) ? 1 : formerValue.getUpperBound();
            if(formerU>bestUpperBoundSoFar) { //avoids the math min, cause we only enter body, if formerValue is too large
                stateUpdate.setQValue(mm.smg.getStatesList().get(s), new StateValue((formerValue==null)?0:formerValue.getLowerBound(), bestUpperBoundSoFar));
                if(setZero){
                    stateUpdate.setZero(mm.smg.getStatesList().get(s),true);
                }
                if (verbose) {
                    mc.getLog().println("Bam! Useful adjustment!");// on state " + smg.getStatesList().get(s) + " former value " + formerValue.getUpperBound() + " new upper bound " + bestUpperBoundSoFar);
                }
            }
        }
    }

    public boolean targetInMec(BitSet mec) throws PrismException{
        for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)) {
            if(stateUpdate.isTarget(mm.smg.getStatesList().get(s))){
                return true;
            }
        }
        return false;
    }

    /**
     * Has its own method, because this is the only difference between this and grey deflate.
     */
    protected double getUpperBoundActionValue(int s, int i) throws PrismException{
        return stateUpdate.getUpperBoundActionValue(mm.smg.getStatesList().get(s),i).getUpperBound();
    }


    /**
     * Update partial model and precompute stuff
     */
    @Override
    public void updatePrecomp() throws PrismException {
        mm.updatePrecomp(seen,mc,stateUpdate);
    }

    /**
     * Stuff only used for children; need to be here, so HeuristicSG can talk about them.
     * Might restructure at some time
     */
    @Override
    public void doPartialBVI() throws PrismException {} //not used for white

    @Override
    public boolean updateCounters(State s, int a, State t) throws PrismException{return false;}

    protected double getLowerBoundActionValue(int s, int i) throws PrismException{
        return stateUpdate.getLowerBoundActionValue(mm.smg.getStatesList().get(s),i).getLowerBound();
    }

}
