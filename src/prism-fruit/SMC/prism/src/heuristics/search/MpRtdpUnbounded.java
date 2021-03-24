//==============================================================================
//
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Mateusz Ujma <mateusz.ujma@cs.ox.ac.uk> (University of Oxford)
//
//------------------------------------------------------------------------------
//
//	This file is part of PRISM.
//
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//==============================================================================

package heuristics.search;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.Distribution;
import explicit.ECComputerFast;
import explicit.MDPSimple;
import explicit.MEC;
import explicit.Model;
import explicit.Product;
import explicit.rewards.MDPRewards;
import heuristics.CachedModelGenerator;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.update.StateUpdate;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2DoubleArrayMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;
import prism.PrismSettings;
import prism.PrismUtils;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.IntConsumer;

// TODO This currently is broken, use CoreLearner
@SuppressWarnings("resource") public class MpRtdpUnbounded
{
	private static final int SAMPLING_VISIT_UPPER_BOUND = 5;
	private static final int REPORT_PROGRESS_EVERY_STEPS = 50000;

	// TODO: Rewrite the whole ModelGenerator thing, properly incorporate state collapsing in state update etc. and separate collapsed model from partial model
	// We have to be careful with step-bounded reachability here - maybe even implement it differently?
	// TODO: MEC detection on collapsed model?

	private final boolean min;
	private final MDPSimple partialMDP = new MDPSimple();
	private final Object2IntMap<State> stateIndexMap = new Object2IntArrayMap<>();
	private final NatBitSet unvisited = NatBitSets.set();
	private final NatBitSet partialModelExploredStates = NatBitSets.set();
	private final State uncertainState;
	private final State minusState;
	private final State plusState;
	private final Map<MEC, MecIterationResult> mecToPreviousSolution = new HashMap<>();
	private final Int2ObjectMap<MEC> stateToMEC = new Int2ObjectArrayMap<>();
	private final NatBitSet viProcessedStates = NatBitSets.set();
	private final Set<State> seen = new HashSet<>();
	private final double meanPayoffPrecision;
	private final double maxRewardUpperBound;
	private final ModelGenerator originalModelGenerator;
	private final MDPRewards partialMDPRewards;
	private final HeuristicsSMGModelChecker smgModelChecker;
	private final CachedModelGenerator cachedModelGenerator;
	private final State initialState;
	private final HeuristicNextState nextState;
	private final StateUpdate stateUpdate;
	private final double tauValue;
	private final ArrayDeque<State> visitedStates = new ArrayDeque<>();
	private final Object2IntArrayMap<State> visitCount = new Object2IntArrayMap<>();
	private int rewardIndex;
	private int bound = -1;
	private boolean verbose;
	//Statistics
	private int trialSteps = 0;
	private int trials = 0;
	private int mecCollapseCount = 0;
	private int mecCollapseThreshold = 1;
	private boolean seenChanged = true;
	private boolean updateSinceLastCollapse = false;
	private int iterationsSinceLastNewState = 0;
	private double[] partialViCurrentResults = new double[0];
	private double[] partialViPreviousResults = new double[0];

	public MpRtdpUnbounded(HeuristicsSMGModelChecker smgModelChecker, StateUpdate stateUpdate, HeuristicNextState nextState,
			CachedModelGenerator cachedModelGenerator, boolean min, double meanPayoffPrecision, double maxRewardUpperBound)
			throws PrismException
	{
		this.smgModelChecker = smgModelChecker;
		this.cachedModelGenerator = cachedModelGenerator;
		this.stateUpdate = stateUpdate;
		this.nextState = nextState;
		this.initialState = cachedModelGenerator.getInitialState();
		this.min = min;

		uncertainState = new State(1);
		uncertainState.setValue(0, "?");

		minusState = new State(1);
		minusState.setValue(0, "-");

		plusState = new State(1);
		plusState.setValue(0, "+");

		partialMDP.setStatesList(new ArrayList<>());
		originalModelGenerator = cachedModelGenerator.getOriginalModelGenerator();
		partialMDPRewards = new PartialMDPRewards();

		this.stateUpdate.setQValue(uncertainState, new StateValue(0d, 1d));
		this.stateUpdate.setQValue(minusState, new StateValue(0d, 0d));
		this.stateUpdate.setQValue(plusState, new StateValue(1d, 1d));
		this.stateUpdate.setTarget(plusState, true);
		this.stateUpdate.setZero(minusState, true);

		this.meanPayoffPrecision = meanPayoffPrecision;
		this.maxRewardUpperBound = maxRewardUpperBound;
		tauValue = smgModelChecker.getSettings().getDouble(PrismSettings.PRISM_MDP_MP_APERIODICITYTAU);
	}

	private double checkReward(double reward, State state) throws PrismException
	{
		if (Double.isNaN(reward)) {
			throw new PrismException(String.format("Reward structure evaluates to NaN at state %s", state));
		}
		if (reward < 0d) {
			throw new PrismException(String.format("Reward structure evaluates to %s at state %s, negative rewards not allowed", reward, state));
		}
		if (reward > maxRewardUpperBound) {
			throw new PrismException(String.format("Reward structure evaluates to %s at state %s, which is bigger than specified upper bound %s",
					reward, state, maxRewardUpperBound));
		}
		return reward;
	}

	/**
	 * The learning algorithm for computing mean-payoff
	 */
	public void computeMeanPayoff() throws PrismException
	{
		State initial = cachedModelGenerator.getInitialState();

		stateUpdate.setQValue(initial, new StateValue(0d, 1d));

		// Start the timer
		long modelCheckingTime = System.currentTimeMillis();

		do {
			exploreAndUpdate(initial);
			trials++;
		} while (!isDone());
		reportProgress();

		// Stop the timer and print log
		long duration = System.currentTimeMillis() - modelCheckingTime;
		smgModelChecker.getLog().println();
		smgModelChecker.getLog().println(String.format("Heuristic model checking time: %.2f secs. Partial model size: %d, %d MECs",
				(double) duration / 1000, partialMDP.getNumStates(), mecToPreviousSolution.size()));
		double lowerBound = getMeanPayoffLowerBound();
		double upperBound = getMeanPayoffUpperBound();
		smgModelChecker.getLog().println(String.format("Result (maximum expected reward): %f, confidence [%f, %f]",
				(lowerBound + upperBound) / 2d, lowerBound, upperBound));
	}

	public double getMeanPayoffLowerBound()
	{
		return getInitialStateValue().lowerBound * maxRewardUpperBound;
	}

	public double getMeanPayoffUpperBound()
	{
		return getInitialStateValue().upperBound * maxRewardUpperBound;
	}

	private void exploreAndUpdate(State state) throws PrismException
	{
		State currentState = state;
		assert visitedStates.isEmpty();
		visitCount.clear();
		visitCount.put(currentState, 0);
		iterationsSinceLastNewState++;

		// EXPLORE Phase
		while (!(currentState.equals(plusState) || currentState.equals(minusState) || currentState.equals(uncertainState))
				&& visitCount.getInt(currentState) < SAMPLING_VISIT_UPPER_BOUND) {
			if (seen.add(currentState)) {
				seenChanged = true;
				iterationsSinceLastNewState = 0;
			}
			visitedStates.addLast(currentState);

			// Update the Q value for the state and fetch the epsilon-close best action
			int bestAction = stateUpdate.update(currentState);
			if (bestAction == -1) {
				smgModelChecker.getLog().println("Ran into a deadlock state " + currentState);
				smgModelChecker.getLog().flush();
				// Deadlock case
				break;
			}
			currentState = nextState.sample(currentState, bestAction);
			visitCount.computeInt(currentState, (mapState, mapIndex) -> mapIndex == null ? 0 : mapIndex + 1);
			trialSteps++;
			if (verbose && trialSteps % REPORT_PROGRESS_EVERY_STEPS == 0) {
				reportProgress();
			}
		}

		if (currentState.equals(uncertainState)) {
			// if (currentState.equals(uncertainState) || currentState.equals(plusState) || currentState.equals(minusState)) { //
			// Sample chose to stay in this MEC. Either it really is good or the bounds are not precise enough. Trigger refinement of the MEC to find out.
			// No need to update the partial model here as a transition can only go to the special states from an already discovered EC. It may be the case
			// that the EC is not maximal, but then exploring the surrounding is more promising than the stay action.
			assert nextState.getPreviousState() != null;
			int stateIndex = stateIndexMap.getInt(nextState.getPreviousState());
			MEC previousStateMEC = stateToMEC.get(stateIndex);
			assert previousStateMEC != null : "Previous state not in a MEC";
			updateMEC(previousStateMEC);
		} else if (!currentState.equals(plusState) && !currentState.equals(minusState)) {
			assert visitedStates.contains(currentState);
			// Iteration stopped because some state was seen multiple times.
			mecCollapseCount++;
			// TODO These values are arbitrary
			if (mecCollapseCount > partialMDP.getNumStates() / 10 && mecCollapseCount > mecCollapseThreshold * 2
					|| iterationsSinceLastNewState > mecCollapseThreshold) {
				// Search for MECs in the partial MDP and process the MECs which were not processed yet
				mecCollapseCount = 0;
				mecCollapseThreshold += 10;
				collapseMECs();
				// Path could now contain multiple states belonging to the same MEC
				visitedStates.clear();
				return;
			}
		} // else: Hit a deadlock

		// UPDATE Phase
		while (!visitedStates.isEmpty()) {
			State visitedState = visitedStates.removeLast();
			// Upper and lower bounds are propagated from the successors
			stateUpdate.update(visitedState);
		}

		trials++;
	}

	public boolean isDone() throws PrismException
	{
		StateValue stateValue = stateUpdate.getQValue(initialState);
		if (stateValue == null) {
			return false;
		}
		double lowerBoundInitialState = stateValue.lowerBound;
		double upperBoundInitialState = stateValue.upperBound;
		return stateUpdate.same(lowerBoundInitialState, upperBoundInitialState);
	}

	private boolean updatePartialModel() throws PrismException
	{
		if (!seenChanged) {
			return false;
		}
		seenChanged = false;

		List<State> statesList = partialMDP.getStatesList();
		IntSet newStates = new IntArraySet();

		// If an initial state doesn't exist in the partial MDP
		// create a new empty state and add it as the initial
		// state. Then map the initial state we know from the
		// cachedModelGenerator, to the index of of the newly added state.
		if (partialMDP.getFirstInitialState() == -1) {
			int index = partialMDP.addState();
			partialMDP.addInitialState(index);
			stateIndexMap.put(initialState, index);
			statesList.add(initialState);
			newStates.add(index);
		}

		for (State seenState : seen) {
			int stateIndex;
			if (!stateIndexMap.containsKey(seenState)) {
				// State is completely new - add it to the partial model
				stateIndex = partialMDP.addState();
				assert !statesList.contains(seenState);
				statesList.add(seenState);
				stateIndexMap.put(seenState, stateIndex);
			} else {
				stateIndex = stateIndexMap.getInt(seenState);
				if (partialModelExploredStates.contains(stateIndex)) {
					// State already explored fully
					assert !unvisited.contains(stateIndex);
					continue;
				} else {
					// This state will be visited later, clear the unvisited flag
					unvisited.clear(stateIndex);
				}
			}
			assert !partialModelExploredStates.contains(stateIndex);
			newStates.add(stateIndex);
		}
		updateSinceLastCollapse |= !newStates.isEmpty();

		if (verbose) {
			smgModelChecker.getLog().println("Update Partial");
			smgModelChecker.getLog().println(String.format("  Updating partial model with %d new states", newStates.size()));
			smgModelChecker.getLog().flush();
		}
		// Iterate over the globally seen states
		IntIterator iterator = newStates.iterator();
		while (iterator.hasNext()) {
			int newStateIndex = iterator.nextInt();
			State newState = partialMDP.getStatesList().get(newStateIndex);
			assert !partialModelExploredStates.contains(newStateIndex);
			partialModelExploredStates.set(newStateIndex);

			List<Distribution> actionsOfStateInPartialMdp = new ArrayList<>();
			originalModelGenerator.exploreState(newState);
			int choices = originalModelGenerator.getNumChoices();
			if (choices == 0) {
				smgModelChecker.getLog().println(String.format("    Warning: Detected a deadlock state %s", newState));
				smgModelChecker.getLog().flush();
				Distribution distribution = new Distribution(newStateIndex, 1.0);
				partialMDP.addChoice(newStateIndex, distribution);
			} else {
				for (int i = 0; i < choices; i++) {
					Distribution distribution = new Distribution();
					int trans = originalModelGenerator.getNumTransitions(i);
					assert trans > 0;
					for (int j = 0; j < trans; j++) {
						double prob = originalModelGenerator.getTransitionProbability(i, j);
						State transitionTarget = originalModelGenerator.computeTransitionTarget(i, j);
						int targetStateIndex = stateIndexMap.computeIntIfAbsent(transitionTarget, state -> {
							int newIndex = partialMDP.addState();
							assert !statesList.contains(state);
							statesList.add(state);
							assert !seen.contains(state) && !unvisited.contains(newIndex);
							unvisited.set(newIndex);
							return newIndex;
						});
						distribution.add(targetStateIndex, prob);
					}
					if (!distribution.isEmpty()) {
						actionsOfStateInPartialMdp.add(distribution);
					}
				}

				for (int actionIndex = 0; actionIndex < actionsOfStateInPartialMdp.size(); actionIndex++) {
					Object action = originalModelGenerator.getChoiceAction(actionIndex);
					Distribution distribution = actionsOfStateInPartialMdp.get(actionIndex);
					int retVal = partialMDP.addActionLabelledChoice(newStateIndex, distribution, action);
					assert retVal != -1;
				}
			}
		}
		// Now: seen = partialModelExploredStates, unvisited are all others
		if (partialViCurrentResults.length < partialMDP.getNumStates()) {
			int newSize = partialMDP.getNumStates() + unvisited.size();
			partialViCurrentResults = Arrays.copyOf(partialViCurrentResults, newSize);
			partialViPreviousResults = Arrays.copyOf(partialViPreviousResults, newSize);
		}

		assert seen.size() + unvisited.size() == partialMDP.getNumStates() :
				String.format("Expected %d states, got %d", seen.size() + unvisited.size(), partialMDP.getNumStates());
		if (verbose) {
			smgModelChecker.getLog().println(String.format("  Number of seen states: %d; Number of unvisited states: %d",
					seen.size(), unvisited.size()));
			smgModelChecker.getLog().flush();
		}
		return true;
	}

	private void collapseMECs() throws PrismException
	{
		updatePartialModel();
		if (!updateSinceLastCollapse) {
			return;
		}
		updateSinceLastCollapse = false;

		long start = System.currentTimeMillis();

		if (verbose) {
			smgModelChecker.getLog().println("Collapse");
			smgModelChecker.getLog().println(String.format("  Searching for MECs in the partial model, %d states", partialModelExploredStates.size()));
			smgModelChecker.getLog().flush();
		}

		// TODO
		List<MEC> mecs = ECComputerFast.computeMECs(partialMDP);

		NatBitSet statesInNewMecs;
		if (mecs.isEmpty()) {
			statesInNewMecs = NatBitSets.set();
		} else {
			statesInNewMecs = NatBitSets.boundedSet(partialMDP.getNumStates());
		}

		List<MEC> newMecs = new ArrayList<>();
		for (MEC mec : mecs) {
			NatBitSet states = mec.states;
			IntIterator iterator = states.iterator();
			while (iterator.hasNext()) {
				int state = iterator.nextInt();
				// Check if any of the states in this MEC have not been touched yet - it is a new MEC then.
				if (!viProcessedStates.contains(state)) {
					statesInNewMecs.or(states);
					newMecs.add(mec);
					break;
				}
			}
		}

		if (verbose) {
			smgModelChecker.getLog().println(String.format("  Found %d MECs, %d new ones, %d states total", mecs.size(), newMecs.size(),
					statesInNewMecs.size()));
			smgModelChecker.getLog().flush();
		}

		// Purge MECs from the solution cache which have expanded
		Iterator<Entry<MEC, MecIterationResult>> solutionCacheIterator = mecToPreviousSolution.entrySet().iterator();
		while (solutionCacheIterator.hasNext()) {
			Entry<MEC, MecIterationResult> solutionCacheEntry = solutionCacheIterator.next();
			NatBitSet mec = solutionCacheEntry.getKey().states;
			if (mec.intersects(statesInNewMecs)) {
				// Some state of this mec is part of the newly explored mecs
				solutionCacheIterator.remove();
				mec.forEach((IntConsumer) state -> {
					partialViCurrentResults[state] = 0d;
					partialViPreviousResults[state] = 0d;
				});
			}
		}

		long duration = System.currentTimeMillis() - start;
		if (verbose) {
			smgModelChecker.getLog().println(String.format("  MEC collapsing done %f secs.", (double) duration / 1000d));
			smgModelChecker.getLog().flush();
		}

		for (MEC mec : newMecs) {
			// Run VI only on the new MECs
			// add the new set of states to CachedModelGenerator
			NatBitSet states = mec.states;
			assert !mecToPreviousSolution.containsKey(mec);
			states.forEach((IntConsumer) mecStateIndex -> stateToMEC.put(mecStateIndex, mec));
			viProcessedStates.or(states);
			updateMEC(mec);
		}
	}

	private void updateMEC(MEC mec) throws PrismException
	{
		MecIterationResult solutionCache = mecToPreviousSolution.computeIfAbsent(mec, bitSet -> new MecIterationResult());
		double requiredPrecision;
		int iterationBound;

		// TODO Adapt the precision and step bound response, add total iteration count to it
		if (solutionCache.getRefinementCount() < 5) {
			// Didn't see this MEC often yet - probably not interesting. Iterate a few times to get a rough idea.
			requiredPrecision = meanPayoffPrecision;
			iterationBound = 10;
		} else {
			// Slowly ramp up the precision and step bound
			requiredPrecision = meanPayoffPrecision / (2d - 1d / (double) (solutionCache.getRefinementCount() - 4));
			iterationBound = 500 * (solutionCache.getRefinementCount() - 4);
		}

		if (solutionCache.getRefinementCount() > 0 && solutionCache.getAchievedPrecision() < requiredPrecision) {
			solutionCache.skipNextIteration();
			// The MEC itself did not change (as a previous solution is present), neither can the actions leaving the MEC change, as all immediate successors
			// where present (potentially as "unvisited" states) in the partial model when the MEC was first discovered.
			return;
		}

		// If MEC has more than 1 state run a step bounded span-seminorm strategy with fixed tau value. Use a increasingly tighter precision to require
		// "arbitrary" precision from the VI. As we only continue the iteration if the ? state is hit and chances of this happening diminish for higher
		// precision, this is fine. Further, don't use adaptive tau, as we aren't allowed to change the MEC dynamically in order to apply the upper and lower
		// bound estimations (at least it is not proven correct).

		double rewardLowerBound;
		double rewardUpperBound;
		Map<State, List<Object2DoubleMap<State>>> actions;
		if (mec.size() == 1) {
			int mecStateIndex = mec.states.firstInt();
			NatBitSet stayingActions = mec.actions.get(mecStateIndex);
			List<Object2DoubleMap<State>> outgoingActions = new ArrayList<>();

			double maximalSelfLoopReward = Double.NEGATIVE_INFINITY;
			for (int choiceIndex = 0; choiceIndex < partialMDP.getNumChoices(mecStateIndex); choiceIndex++) {
				if (stayingActions.contains(choiceIndex)) {
					assert partialMDP.allSuccessorsInSet(mecStateIndex, choiceIndex, mec.states);
					maximalSelfLoopReward = Math.max(maximalSelfLoopReward, partialMDPRewards.getTransitionReward(mecStateIndex, choiceIndex));
				} else {
					assert !partialMDP.allSuccessorsInSet(mecStateIndex, choiceIndex, mec.states);
					Object2DoubleMap<State> distribution = getDistributionInPartialModel(mecStateIndex, choiceIndex);
					assert !distribution.containsKey(plusState)
							|| Arrays.asList(plusState, minusState, uncertainState).containsAll(distribution.keySet())
							: String.format("State %d has a mixed transition (%s)", mecStateIndex, distribution);
					if (!distribution.isEmpty()) {
						outgoingActions.add(distribution);
					}
				}
			}

			if (stayingActions.isEmpty()) {
				rewardLowerBound = 0d;
				rewardUpperBound = 0d;
			} else {
				assert maximalSelfLoopReward >= 0;
				maximalSelfLoopReward += partialMDPRewards.getStateReward(mecStateIndex);
				rewardLowerBound = maximalSelfLoopReward;
				rewardUpperBound = maximalSelfLoopReward;
			}
			solutionCache.setNextIteration(rewardUpperBound, rewardLowerBound, requiredPrecision);

			if (!outgoingActions.isEmpty()) {
				State mecState = partialMDP.getStatesList().get(mecStateIndex);
				actions = Collections.singletonMap(mecState, outgoingActions);
			} else {
				actions = Collections.emptyMap();
			}
		} else {
			PartialViResults partialViResults = runVISteps(mec, requiredPrecision, iterationBound);
			rewardLowerBound = partialViResults.getLowerBound();
			rewardUpperBound = partialViResults.getUpperBound();
			solutionCache.setNextIteration(rewardLowerBound, rewardUpperBound, requiredPrecision);

			// Collect all the actions leaving the MEC
			actions = getAllLeavingMECinPartial(mec);

			if (solutionCache.getRefinementCount() > 2 && verbose) {
				smgModelChecker.getLog().println(
						String.format("New bounds for non-trivial MEC: [%.15f, %.15f]. %d states with outgoing edges (%d transitions total)",
								rewardLowerBound, rewardUpperBound, actions.size(), actions.values().stream().mapToInt(List::size).sum()));
				smgModelChecker.getLog().flush();
			}
		}

		List<Object2DoubleMap<State>> allActions = new ArrayList<>();
		actions.values().forEach(allActions::addAll);

		// Add edges to +, - and ? states
		// Dividing by the maximal reward so that we have proper probabilities
		Object2DoubleMap<State> stayAction = new Object2DoubleArrayMap<>();
		assert rewardLowerBound / maxRewardUpperBound <= 1d;
		assert (1d - (rewardUpperBound / maxRewardUpperBound)) <= 1d;
		if (rewardLowerBound > 0d) {
			stayAction.put(plusState, rewardLowerBound / maxRewardUpperBound);
		}
		if (rewardUpperBound < maxRewardUpperBound) {
			stayAction.put(minusState, 1d - rewardUpperBound / maxRewardUpperBound);
		}
		if (rewardLowerBound < rewardUpperBound) {
			stayAction.put(uncertainState, (rewardUpperBound - rewardLowerBound) / maxRewardUpperBound);
		}
		allActions.add(stayAction);

		Set<State> states = new HashSet<>();

		mec.states.forEach((IntConsumer) j -> {
			State state = partialMDP.getStatesList().get(j);
			assert !Arrays.asList(plusState, minusState, uncertainState).contains(state);
			states.add(state);
		});

		cachedModelGenerator.mergeStates(states, allActions);
	}

	private PartialViResults runVISteps(MEC mec, double precision, int maxSteps)
	{
		double[] previousValues = partialViPreviousResults;
		double[] currentValues = partialViCurrentResults;

		NatBitSet states = mec.states;
		int firstState = states.firstInt();
		// Start iterations
		int iteration = 0;
		double maximum = Double.NaN;
		double minimum = Double.NaN;
		while (iteration < maxSteps) {
			double[] finalPreviousValues = previousValues;
			partialMDP.mvMultRewMinMaxAperiodic(s -> finalPreviousValues[s], partialMDPRewards, min, currentValues, mec, null, tauValue);
			iteration++;

			maximum = 0d;
			minimum = Double.POSITIVE_INFINITY;

			IntIterator iterator = states.iterator();
			while (iterator.hasNext()) {
				int state = iterator.nextInt();
				double difference = currentValues[state] - previousValues[state];
				assert difference >= 0d || PrismUtils.doublesAreClose(difference, 0d, 1e-8, true) :
						String.format("Negative difference %.10f in state %d, first %d", difference, state, firstState);
				maximum = Math.max(maximum, difference);
				minimum = Math.min(minimum, difference);
			}
			assert !Double.isNaN(maximum) && !Double.isNaN(minimum);
			double norm = maximum - minimum;

			assert norm >= 0d;
			if (norm < tauValue * precision) {
				break;
			}

			if (currentValues[firstState] > 1000d) {
				double referenceValue = previousValues[firstState];
				IntIterator stateIterator = states.iterator();
				while (stateIterator.hasNext()) {
					int state = stateIterator.nextInt();
					currentValues[state] = currentValues[state] - referenceValue;
					previousValues[state] = previousValues[state] - referenceValue;
				}
			}
			double[] swap = previousValues;
			previousValues = currentValues;
			currentValues = swap;
		}
		//noinspection ArrayEquality This is intended
		if (currentValues != partialViCurrentResults) {
			partialMDP.mvMultRewMinMaxAperiodic(s -> partialViPreviousResults[s], partialMDPRewards, min, partialViCurrentResults, mec, null, tauValue);
		}

		assert !Double.isNaN(maximum) && !Double.isNaN(minimum) && minimum <= maximum;
		return new PartialViResults(minimum / tauValue, maximum / tauValue);
	}

	private Map<State, List<Object2DoubleMap<State>>> getAllLeavingMECinPartial(MEC mec)
	{
		Map<State, List<Object2DoubleMap<State>>> actions = new HashMap<>();
		List<State> stateList = partialMDP.getStatesList();
		NatBitSet states = mec.states;

		IntIterator stateIterator = states.iterator();
		while (stateIterator.hasNext())  {
			int stateInMec = stateIterator.nextInt();
			// Go through each state of the MEC and collect its distributions
			List<Object2DoubleMap<State>> distributions = new ArrayList<>();
			for (int action = 0; action < partialMDP.getNumChoices(stateInMec); action++) {
				// Add the distribution only if not all of the successors are within the MEC
				double outgoingWeight = 0d;
				Object2DoubleMap<State> distribution = new Object2DoubleArrayMap<>();
				for (Int2DoubleMap.Entry entry : partialMDP.getTransitions(stateInMec, action)) {
					if (states.contains(entry.getIntKey())) {
						continue;
					}
					outgoingWeight += entry.getDoubleValue();
					distribution.put(stateList.get(entry.getIntKey()), entry.getDoubleValue());
				}
				if (outgoingWeight == 0d) {
					continue;
				}
				Object2DoubleMap<State> weightedDistribution;
				if (outgoingWeight == 1d) {
					weightedDistribution = distribution;
				} else {
					weightedDistribution = new Object2DoubleArrayMap<>(distribution.size());
					for (Object2DoubleMap.Entry<State> entry : distribution.object2DoubleEntrySet()) {
						weightedDistribution.put(entry.getKey(), entry.getDoubleValue() / outgoingWeight);
					}
				}
				assert !weightedDistribution.containsKey(plusState)
						|| Arrays.asList(plusState, minusState, uncertainState).containsAll(weightedDistribution.keySet())
						: "State " + stateInMec + " has a mixed transition (" + weightedDistribution + ")";
				if (weightedDistribution.containsKey(plusState)) {
					continue;
				}
				distributions.add(weightedDistribution);
			}
			if (distributions.isEmpty()) {
				continue;
			}
			State sourceState = stateList.get(stateInMec);
			actions.put(sourceState, distributions);
		}
		return actions;
	}

	private Object2DoubleMap<State> getDistributionInPartialModel(int stateIndex, int choiceIndex)
	{
		Object2DoubleMap<State> distribution = new Object2DoubleArrayMap<>();

		for (Int2DoubleMap.Entry transitionEntry : partialMDP.getTransitions(stateIndex, choiceIndex)) {
			distribution.put(partialMDP.getStatesList().get(transitionEntry.getIntKey()), transitionEntry.getDoubleValue());
		}

		return distribution;
	}

	public void setRewardIndex(int rewardIndex)
	{
		this.rewardIndex = rewardIndex;
	}

	public int getBound()
	{
		return bound;
	}

	public void setBound(int bound)
	{
		this.bound = bound;
	}

	public void setVerbose(boolean v)
	{
		this.verbose = v;
	}

	public boolean getVerbose()
	{
		return this.verbose;
	}

	public StateValue getInitialStateValue()
	{
		StateValue qValue = stateUpdate.getQValue(initialState);
		if (qValue == null) {
			return new StateValue(-1.0, -1.0);
		}
		return qValue;
	}

	private void reportProgress()
	{
		if (verbose) {
			String progressString = String.format("%n== Progress Report ==%n"
							+ "  Trials: %d, steps: %d, avg len: %f%n"
							+ "  States: %d in partial model, %d with value, %d unvisited%n"
							+ "  Bounds [%.15f,%.15f], diff: %.6g%n",
					trials, trialSteps, (double) trialSteps / (double) trials,
					partialMDP.getNumStates(), stateUpdate.getNumExplored(), unvisited.size(),
					getMeanPayoffLowerBound(), getMeanPayoffUpperBound(), getMeanPayoffUpperBound() - getMeanPayoffLowerBound());
			smgModelChecker.getLog().println(progressString);
			smgModelChecker.getLog().flush();
		}
	}

	private static final class PartialViResults
	{
		private final double lowerBound;
		private final double upperBound;

		private PartialViResults(double lowerBound, double upperBound)
		{
			this.lowerBound = lowerBound;
			this.upperBound = upperBound;
		}

		public double getLowerBound()
		{
			return lowerBound;
		}

		public double getUpperBound()
		{
			return upperBound;
		}
	}

	private static class MecIterationResult
	{
		private int refinementCount = 0;
		private double upperBound = -1d;
		private double lowerBound = -1d;
		private double precision;

		void setNextIteration(double lowerBound, double upperBound, double precision)
		{
			assert 0d <= lowerBound && lowerBound <= upperBound : String.format("Invalid bounds [%f,%f]", lowerBound, upperBound);
			this.precision = precision;
			refinementCount += 1;
			this.upperBound = upperBound;
			this.lowerBound = lowerBound;
		}

		double getPrecision()
		{
			return precision;
		}

		int getRefinementCount()
		{
			return refinementCount;
		}

		double getAchievedPrecision()
		{
			assert refinementCount > 0;
			assert 0d <= lowerBound && lowerBound <= upperBound;
			return upperBound - lowerBound;
		}

		void skipNextIteration()
		{
			assert refinementCount > 0;
			refinementCount++;
		}
	}

	private class PartialMDPRewards implements MDPRewards
	{
		private double[] stateRewardCache = new double[0];
		private double[][] transitionRewardCache = new double[0][];

		private double getStateRewardFromCache(int stateIndex)
		{
			if (stateRewardCache.length <= stateIndex) {
				int previousSize = stateRewardCache.length;
				stateRewardCache = Arrays.copyOf(stateRewardCache, partialMDP.getStatesList().size() * 2);
				Arrays.fill(stateRewardCache, previousSize, stateRewardCache.length, Double.NaN);
				double reward = getStateRewardFromModel(stateIndex);
				stateRewardCache[stateIndex] = reward;
				return reward;
			}
			double cachedReward = stateRewardCache[stateIndex];
			if (Double.isNaN(cachedReward)) {
				double reward = getStateRewardFromModel(stateIndex);
				stateRewardCache[stateIndex] = reward;
				return reward;
			}
			return cachedReward;
		}

		private double getStateRewardFromModel(int stateIndex)
		{
			State state = partialMDP.getStatesList().get(stateIndex);
			try {
				return checkReward(originalModelGenerator.getStateReward(rewardIndex, state), state);
			} catch (PrismException e) {
				// PrismException not allowed here
				throw new AssertionError(e);
			}
		}

		private double getTransitionRewardFromModel(int stateIndex, int actionIndex)
		{
			State state = partialMDP.getStatesList().get(stateIndex);
			Object action = partialMDP.getAction(stateIndex, actionIndex);
			try {
				return checkReward(originalModelGenerator.getStateActionReward(rewardIndex, state, action), state);
			} catch (PrismException e) {
				// PrismException not allowed here
				throw new AssertionError(e);
			}
		}

		private double getTransitionRewardFromCache(int stateIndex, int transitionIndex)
		{
			if (transitionRewardCache.length <= stateIndex) {
				transitionRewardCache = Arrays.copyOf(transitionRewardCache, partialMDP.getStatesList().size() * 2);
			}
			double[] transitionArray = transitionRewardCache[stateIndex];
			if (transitionArray == null) {
				double[] newTransitionArray = new double[partialMDP.getNumChoices(stateIndex)];
				Arrays.fill(newTransitionArray, Double.NaN);
				transitionRewardCache[stateIndex] = newTransitionArray;
				double reward = getTransitionRewardFromModel(stateIndex, transitionIndex);
				newTransitionArray[transitionIndex] = reward;
				return reward;
			}
			double cachedReward = transitionArray[transitionIndex];
			if (Double.isNaN(cachedReward)) {
				double reward = getTransitionRewardFromModel(stateIndex, transitionIndex);
				transitionArray[transitionIndex] = reward;
				return reward;
			}
			return cachedReward;
		}

		@Override public double getStateReward(int stateIndex)
		{
			assert partialModelExploredStates.contains(stateIndex);
			double reward = getStateRewardFromCache(stateIndex);
			assert reward == getStateRewardFromModel(stateIndex);
			return reward;
		}

		@Override public double getTransitionReward(int stateIndex, int actionIndex)
		{
			assert partialModelExploredStates.contains(stateIndex);
			double reward = getTransitionRewardFromCache(stateIndex, actionIndex);
			assert reward == getTransitionRewardFromModel(stateIndex, actionIndex);
			return reward;
		}

		@Override public MDPRewards liftFromModel(Product<? extends Model> product)
		{
			throw new UnsupportedOperationException("Lift not supported");
		}

		@Override public boolean hasTransitionRewards()
		{
			return true;
		}
	}
}