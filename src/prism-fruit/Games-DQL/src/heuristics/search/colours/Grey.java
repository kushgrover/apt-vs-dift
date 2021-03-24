package heuristics.search.colours;

import explicit.ModelExplorer;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.search.StateValue;
import heuristics.update.StateUpdate;
import parser.State;
import prism.Pair;
import prism.PrismException;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

public class Grey extends White {
    protected int statesOA;
    protected int availableOA;
    protected double delta;
    protected double pmin;
    protected int postMax;

    protected double transDelta;

    public Grey(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min, int statesOA, int availableOA, double delta, double pmin, int postMax) throws PrismException {
        super(mc, su, nextState, pme, min);
        this.statesOA=statesOA;
        this.availableOA=availableOA;
        this.delta=delta;
        this.pmin=pmin;
        this.postMax=postMax;

        mc.getLog().println("IsMDP " + mm.isMDP + " Collapse " + collapseMECs + " break " + breakImmediately);
        mc.getLog().println("ColourParams: S:" + statesOA + " Av: " + availableOA + " eps: " +stateUpdate.getEpsilon() + " delta: " + delta + " pmin: " + pmin);
        //We distribute the allowed error over all transitions; a priori, we can only use the states and actions and postmax overapproximations
        //For partialBVI, we can do better (see method doPartialBVI)
        //Since in phase we don't give any guarantees at all, we could also just use delta or whatever we want; let's see what behaves best
        this.transDelta = delta/(statesOA * availableOA * postMax);
        mc.getLog().println("TransDelta: " + transDelta);
        specialStuff(); //needed to differentiate between constructor of grey and black
    }
    public void specialStuff(){
        this.oc = new OccurrenceCounter(availableOA, transDelta, pme);
    }

    protected static final Random random = new Random();
    @Override
    public int sampleAction(State currentState) throws PrismException{
        ArrayList<Integer> bestActions = getBestActions(currentState);
        //pick a random action from the set of best actions

        int bestAction = bestActions.get(random.nextInt(bestActions.size()));
        return bestAction;
    }

    public ArrayList<Integer> getBestActions(State s) throws PrismException{
        pme.queryState(s);
        int choices = pme.getNumChoices();//ok even in black, since we have oracle for actions
        ArrayList<Integer> bestActions = new ArrayList<>();
        if(choices == 0){
            bestActions.add(-1);
        }
        if(choices==1){//If there is only one action, use it without computing any more stuff
            bestActions.add(0);
            return bestActions;
        }
        boolean isMaxState = mm.isMaxState(s);
        //Random
//        for (int i = 0; i < choices; i++){
//            bestActions.add(i);
//        }

//        //Sensible thing
        double currentBound;
        double bestBound = isMaxState ? 0 : 1;
        for(int a=0;a<choices;a++) {
            currentBound = isMaxState ? getUpperBoundActionValue(s,a) : getLowerBoundActionValue(s, a);
            //Also take into account the confidence width TODO: Is this reasonable and helping us?
            double confWidth = oc.getConfidenceForSApair(s,a) * 10; //times 10, to make confWidth more important
            if(isMaxState){
                //Large confWidth => Larger value => More attractive
                currentBound = currentBound * confWidth;
            }
            else{
                //Large confWidth => Small inverse => Small value => More attractive
                currentBound = Math.min(currentBound * (1.0/confWidth),1);
            }
            if(isBetter(currentBound,bestBound,isMaxState)){
                //new best Action, remove all others we had until now
                bestActions.clear();
                bestActions.add(a);
                bestBound = currentBound;
            }
            else if (!isBetter(bestBound,currentBound,isMaxState)){
                //Action is equally good, add it to set
                bestActions.add(a);
            }
        }
        if(bestActions.isEmpty()){throw new PrismException("Empty set of best actions in Grey. This should never happen.");}

        return bestActions;
    }

    @Override
    public boolean updateCounters(State s, int a, State t){
        oc.incrementTriple(s,a,t);
        return false;
    }


    @Override
    public void update() throws PrismException{
        State lastState = visited.removeLast();//remove last state; if it is a target, we do not have any probabilities for it; otherwise we are in a loop and it is still in the path somewhere else
        if (stateUpdate.getQValue(lastState)==null && stateUpdate.isTarget(lastState)){
            stateUpdate.setQValue(lastState, new StateValue(1,1)); //need to set the value for each target once
        }
        oc.updateAllProbs(new HashSet(visited));
        while (!visited.isEmpty()) {
            State visitedState = visited.removeLast();
            /** Update of QValue implemented here in Grey, because otherwise I have to update both StateUpdateMDP and -SG, and I need more methods and pass things around...
             * Hacky, but easier.*/
            StateValue currVal = getBestValue(visitedState);
            StateValue formerVal = stateUpdate.getQValue(visitedState);
            if(formerVal==null){
                stateUpdate.setQValue(visitedState,currVal);
            }
            else{//make monotonic
                stateUpdate.setQValue(visitedState,new StateValue(Math.max(currVal.getLowerBound(),formerVal.getLowerBound()),Math.min(currVal.getUpperBound(),formerVal.getUpperBound())));
            }
        }
    }

    /**
     * @param s State we want to have the bound for
     * @return The best bound we can estimate at the moment, using lowerProbEstimates and remainingProb from oc
     */
    protected StateValue getBestValue(State s) throws PrismException{
        if (stateUpdate.isZero(s)){
            return stateUpdate.getQValue(s);
        }
        double bestLower, bestUpper; //initialized to the worst possible value for the player
        boolean max = mm.isMaxState(s); //used for initializing current best bounds as well as to compare whether sth is better
        if(max) {
            bestLower = 0;
            bestUpper = 0;
        }
        else{
            bestLower = 1;
            bestUpper = 1;
        }

        //It is ok to use pme even for black, since we assume to have an oracle for available actions
        pme.queryState(s);
        int choices=pme.getNumChoices();
        if(choices==0){
            return stateUpdate.isTarget(s) ? new StateValue(1,1) : new StateValue(0,0);
        }

        for (int a=0; a<choices;a++){
            Pair<Double,Double> bounds = getEstimatedActionBound(s, a);
            double Lsa = bounds.first;
            double Usa = bounds.second;
            if (isBetter(Lsa, bestLower, max)){
                bestLower = Lsa;
            }
            if (isBetter(Usa, bestUpper, max)){
                bestUpper = Usa;
            }
        }
        return new StateValue(bestLower,bestUpper);
    }

    /**
     * Convenience to make things more readable.
     * Returns whether b1 is better than b2 for a player; if max is true, larger is better, else smaller.
     */
    protected boolean isBetter(double b1, double b2, boolean max){
        if(max){
            return b1 > b2;
        }
        else{
            return b1 < b2;
        }
    }

    protected Pair<Double,Double> getEstimatedActionBound(State s, int a) throws PrismException{
        //handle base cases
        if(!oc.allSeenStates().contains(s) || oc.getPairCount(s,a)==0){
            //First time we see the state; or an action that has not been played yet. We know nothing.
            return new Pair<>(0d,1d);
        }

        double lowerEst = 0;
        double upperEst = 0;
        //Next three vars for the conservative guess in the end.
        double[] allSuccBoundsL = new double[oc.getSeenSuccs(s,a).size()];
        double[] allSuccBoundsU = new double[oc.getSeenSuccs(s,a).size()];
        int i = 0;
        for (State t : oc.getSeenSuccs(s,a)){
            StateValue val = stateUpdate.getQValue(t);
            double tUBound, tLBound;
            //handles that val could be null
            if(val==null){
                tLBound = 0;
                tUBound = 1;
            }
            else{
                tLBound = val.getLowerBound();
                tUBound = val.getUpperBound();

            }
            allSuccBoundsL[i] = tLBound;//remember bounds for conservative guess
            allSuccBoundsU[i] = tUBound;//remember bounds for conservative guess
            i++;
            double probEst = oc.getLowerProbEstimate(s,a,t);
            lowerEst += tLBound * probEst;
            upperEst += tUBound * probEst;
        }
        pme.queryState(s);
        double remProb = oc.getRemainingProb(s,a);
        lowerEst +=  remProb * conservativeGuess(true,allSuccBoundsL,pme.getNumChoices());
        upperEst +=  remProb * conservativeGuess(false,allSuccBoundsU,pme.getNumChoices());
        return new Pair<>(lowerEst,upperEst);
    }

    protected double conservativeGuess(boolean lower, double[] allSuccBounds, int numSuccs) throws PrismException{
        if(allSuccBounds.length< numSuccs){
            return lower ? 0 : 1; //gotta catch this, because if we have not seen all succs (or none at all) we return bullshit otherwise
        }

        double mostConservative;
        if (lower){
            mostConservative = 1; //worst possible value
        }
        else{
            mostConservative = 0; //worst possible value
        }
        for (int i=0; i<allSuccBounds.length; i++){
            mostConservative = isBetter(allSuccBounds[i], mostConservative, !lower) ? allSuccBounds[i] : mostConservative;
        }
        return mostConservative;
    }


    //Handling BECs and precomp stays the same, except that we need to use the prob estimate to get best leaving U
    @Override
    protected double getUpperBoundActionValue(int s, int i) throws PrismException{
        return getEstimatedActionBound(mm.smg.getStatesList().get(s), i).second;
    }

    //This version is used for finding best action in explore
    protected double getUpperBoundActionValue(State s, int i) throws PrismException{
        return getEstimatedActionBound(s, i).second;
    }

    protected double getLowerBoundActionValue(State s, int i) throws PrismException{
        return getEstimatedActionBound(s, i).first;
    }


    @Override
    public void doPartialBVI() throws PrismException{
        oc.updateAllProbs(oc.allSeenStates());
        //mc.getLog().println("Running precomputation...");
        updatePrecomp();

        //we need to give some error probability to those transitions where we could fail.
        //so we count all transitions in the partial model, excluding actions with only one successor, since we can't fail there (in grey)
        int numTrans = 0;
        //get transDelta for current partial model
        for (int s = 0; s<mm.smg.getStatesList().size(); s++){
            for (int a = 0; a<mm.smg.getNumChoices(s); a++){
                int trans = mm.smg.getNumTransitions(s,a);
                if(trans==1){
                    continue;
                }
                else{
                    numTrans += trans;
                }
            }
        }
        oc.setTransDelta(delta/numTrans);

        PartialBVI.progressEpsilon=stateUpdate.getEpsilon();

        //Note: unvisited contains the rim states after we buildPartialModel
        PartialBVI.partialBVI(oc,mm.unvisited,stateUpdate,mm.smg,min,mc.getLog(), false);
        printProgress(System.currentTimeMillis() - modelCheckingTime);

        oc.setTransDelta(transDelta); //set back to old value, whatever that was before
    }
}
