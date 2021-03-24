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
import explicit.ConstructModel;
import explicit.Distribution;
import explicit.MDP;
import explicit.MDPModelChecker;
import explicit.MDPSimple;
import heuristics.CachedModelGenerator;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.update.StateUpdate;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.objects.Object2DoubleArrayMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import parser.State;
import parser.ast.Expression;
import prism.ModelGenerator;
import prism.PrismException;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.IntConsumer;

public class HeuristicRtdpUnbounded extends Heuristic
{
	private static final int PROGRESS_REPORT_TRIALS = 1000;

	private final MDPSimple partialMDP = new MDPSimple();
	private final Object2IntMap<State> stateIndexMap = new Object2IntOpenHashMap<>();
	private final NatBitSet unvisited = NatBitSets.set();
	private final NatBitSet partialModelExploredStates = NatBitSets.set();
	private final Set<State> seen = new HashSet<>();
	private final Deque<State> visited = new ArrayDeque<>();
	private long modelCheckingTime = 0L;
	private final int mecCollapseThreshold = 10;
	private int mecCollapseCount = 0;
	private boolean seenChanged = false;
	private boolean updateSinceLastCollapse = false;
	private final NatBitSet mecStates = NatBitSets.set();

	public HeuristicRtdpUnbounded(HeuristicsSMGModelChecker smgModelChecker, StateUpdate stateUpdate, HeuristicNextState nextState,
			ModelGenerator modelGenerator, boolean min)
			throws PrismException
	{
		super(smgModelChecker, stateUpdate, nextState, modelGenerator, min);
	}

	@Override protected void heuristicStart() throws PrismException
	{
		modelCheckingTime = System.currentTimeMillis();
	}

	@Override protected void heuristicStop() throws PrismException
	{
		long duration = System.currentTimeMillis() - modelCheckingTime;
		smgModelChecker.getLog().println();
		smgModelChecker.getLog().println("Heuristic model checking time in " + ((double) duration / 1000) + " secs.");
	}

	@Override public void exploreAndUpdate(State state) throws PrismException
	{
		State currentState = state;

		// EXPLORE Phase
		while (!stateUpdate.isTarget(currentState) && !stateUpdate.isZeroState(currentState) && !visited.contains(currentState)) {
			// Update the Q value for the state and fetch the epsilon-close best action
			seenChanged |= seen.add(currentState);
			visited.addLast(currentState);

			int bestAction = stateUpdate.update(currentState);
			currentState = nextState.sample(currentState, bestAction);

			// precomp code was here, see e.g. commit dc44b4f01e7cb46fb249993186bf55377da530d7
			trialSteps++;
		}

		stateUpdate.setTargetOrZeroValue(currentState);

		if (visited.contains(currentState)) {
			mecCollapseCount++;
			if (mecCollapseCount > mecCollapseThreshold) {
				// Construct partialMDP from all the states seen so far
				updatePrecomp();

				// Search for MECs in the partial MDP and collapse them
				collapseMECs();
				mecCollapseCount = 0;
			}
		}

		// UPDATE Phase
		while (!visited.isEmpty()) {
			State visitedState = visited.removeLast();
			// Upper and lower bounds are propagated from the successors
			stateUpdate.update(visitedState);
		}

		if ((trials + 1) % PROGRESS_REPORT_TRIALS == 0) {
			reportProgress(trials, trialSteps);
		}
		trials++;
	}

	@Override public boolean isDone() throws PrismException
	{
		StateValue stateValue = stateUpdate.getQValue(initialState);
		if (stateValue != null) {
			double lowerBoundInitialState = stateValue.lowerBound;
			double upperBoundInitialState = stateValue.upperBound;
			return stateUpdate.same(lowerBoundInitialState, upperBoundInitialState);
		}
		return false;
	}

	private void updatePrecomp() throws PrismException
	{
		buildPartialModel();
		partialMDP.findDeadlocks(true);

		MDPModelChecker mdpmc = getMC();
		NatBitSet target = computeTarget(mdpmc, partialMDP);

		NatBitSet t = NatBitSets.set();
		t.or(target);
		t.or(unvisited);

		if (min) {
			// Compute states which have a 0 probability of reaching the target
			NatBitSet prob0 = computeProb0(mdpmc, partialMDP, t);

			prob0.forEach((IntConsumer) i -> {
				State s = partialMDP.getStatesList().get(i);
				stateUpdate.setQValue(s, new StateValue(0, 0));
				stateUpdate.setZero(s, true);
			});

			// Compute states which can reach the target with probability 1
			NatBitSet prob1 = computeProb1(mdpmc, partialMDP, target);

			prob1.forEach((IntConsumer) i -> {
				State s = partialMDP.getStatesList().get(i);
				if (seen.contains(s)) {
					stateUpdate.setQValue(s, new StateValue(1, 1));
					stateUpdate.setTarget(s, true);
				}
			});
		} else {
			// Compute states which have a 0 probability of reaching the target
			NatBitSet prob0 = computeProb0(mdpmc, partialMDP, t);

			prob0.forEach((IntConsumer) i -> {
				State s = partialMDP.getStatesList().get(i);
				if (seen.contains(s)) {
					stateUpdate.setQValue(s, new StateValue(0, 0));
					stateUpdate.setZero(s, true);
				}
			});
		}
	}

	private void buildPartialModel() throws PrismException
	{
		if (!seenChanged) {
			return;
		}
		seenChanged = false;

		Set<State> newStates = new HashSet<>();

		List<State> statesList = partialMDP.getStatesList();
		if (statesList == null) {
			statesList = new ArrayList<>();
			partialMDP.setStatesList(statesList);
		}

		unvisited.clear();

		if (partialMDP.getFirstInitialState() == -1) {
			int index = partialMDP.addState();
			partialMDP.addInitialState(index);
			stateIndexMap.put(initialState, index);
			statesList.add(initialState);
			newStates.add(initialState);
		}

		for (State seenState : seen) {
			if (!stateIndexMap.containsKey(seenState)) {
				// State is completely new - add it to the partial model
				int newStateIndex = partialMDP.addState();
				assert !statesList.contains(seenState);
				statesList.add(seenState);
				stateIndexMap.put(seenState, newStateIndex);
			} else {
				int stateIndex = stateIndexMap.getInt(seenState);
				if (partialModelExploredStates.contains(stateIndex)) {
					// State already explored fully
					assert !unvisited.contains(stateIndex);
					continue;
				}
				// This state will be visited later, clear the unvisited flag
				assert !partialModelExploredStates.contains(stateIndex);
				unvisited.clear(stateIndex);
			}
			newStates.add(seenState);
		}
		updateSinceLastCollapse |= !newStates.isEmpty();

		if (verbose) {
			smgModelChecker.getLog().println(String.format("Updating partial model with %d new states", newStates.size()));
			smgModelChecker.getLog().flush();
		}
		// Iterate over the globally seen states
		for (State newState : newStates) {
			int index = stateIndexMap.getInt(newState);
			assert !partialModelExploredStates.contains(index);
			partialModelExploredStates.set(index);
			List<Distribution> dists = buildAllDistributionsIn(partialMDP, newState, seen, statesList);
			partialMDP.clearState(index);
			for (Distribution dist : dists) {
				partialMDP.addChoice(index, dist);
			}
		}

		smgModelChecker.getLog().println("Model built " + partialMDP.getNumStates());
		smgModelChecker.getLog().flush();
	}

	private List<Distribution> buildAllDistributionsIn(MDPSimple mdp, State s, Set<State> seen, List<State> statesList) throws PrismException
	{
		List<Distribution> dists = new ArrayList<>();
		modelGenerator.exploreState(s);
		int choices = modelGenerator.getNumChoices();
		for (int i = 0; i < choices; i++) {
			Distribution d = new Distribution();
			int trans = modelGenerator.getNumTransitions(i);
			boolean inSeen = false;
			for (int j = 0; j < trans; j++) {
				double prob = modelGenerator.getTransitionProbability(i, j);
				State t = modelGenerator.computeTransitionTarget(i, j);
				inSeen = inSeen || seen.contains(t);
				int index;
				if (stateIndexMap.containsKey(t)) {
					index = stateIndexMap.getInt(t);
				} else {
					index = mdp.addState();
					stateIndexMap.put(t, index);
					statesList.add(t);
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

	private void preComp()
	{
		try {
			ConstructModel constructModel = new ConstructModel(smgModelChecker);
			// TODO: PORT, this is a place where things can go wrong
			MDP mdp = (MDP) constructModel.constructModel(modelGenerator);
			MDPModelChecker mdpModelChecker = getMC();
			NatBitSet target = computeTarget(mdpModelChecker, mdp);
			NatBitSet zero = computeProb0(mdpModelChecker, mdp, target);
			NatBitSet one = computeProb1(mdpModelChecker, mdp, target);

			List<State> statesList = mdp.getStatesList();
			for (int i = 0; i < statesList.size(); i++) {
				State s = statesList.get(i);
				if (target.contains(i) || one.contains(i)) {
					stateUpdate.setTarget(s, true);
				} else {
					stateUpdate.setTarget(s, false);
				}
				if (zero.contains(i)) {
					stateUpdate.setZero(s, true);
				}
				modelGenerator.exploreState(s);
			}
			if (!min) {
				collapseMECs();
			}
		} catch (PrismException e) {
			smgModelChecker.getLog().print(e);
		}
	}

	private NatBitSet computeProb0(MDPModelChecker mdpModelChecker, MDP mdp, NatBitSet target)
	{
		return mdpModelChecker.prob0(mdp, null, target, min, null);
	}

	private NatBitSet computeProb1(MDPModelChecker mdpModelChecker, MDP mdp, NatBitSet target)
	{
		return mdpModelChecker.prob1(mdp, null, target, min, null);
	}

	private NatBitSet computeTarget(MDPModelChecker mdpModelChecker, MDP mdp) throws PrismException
	{
		return mdpModelChecker.checkExpression(mdp, stateUpdate.getTarget(), null).getNatBitSet();
	}

	private MDPModelChecker getMC() throws PrismException
	{
		MDPModelChecker mdpModelChecker = new MDPModelChecker(smgModelChecker);
		mdpModelChecker.setModulesFileAndPropertiesFile(smgModelChecker.getModulesFile(), smgModelChecker.getPropertiesFile(), null);
		return mdpModelChecker;
	}

	private void collapseMECs() throws PrismException
	{
		if (!updateSinceLastCollapse) {
			return;
		}

		smgModelChecker.getLog().println("Starting MEC collapsing");
		smgModelChecker.getLog().flush();
		long start = System.currentTimeMillis();
		NatBitSet all = NatBitSets.set();
		NatBitSet notOne = NatBitSets.set();
		for (int i = 0; i < partialMDP.getNumStates(); i++) {
			all.set(i);
			StateValue sv = stateUpdate.getQValue(partialMDP.getStatesList().get(i));
			if (sv != null) {
				if (sv.upperBound < 1.0) {
					notOne.set(i);
				}
			}
		}
		all.xor(notOne);
		// Don't bother with unvisited states - if called with no argument, the method allocates a bit set anyway
		all.and(partialModelExploredStates);
		explicit.ECComputer ec = explicit.ECComputer.createECComputer(smgModelChecker, partialMDP);
		ec.computeMECStates(all);
		List<NatBitSet> mecs = ec.getMECStates();

		List<NatBitSet> newMecs = new ArrayList<>();
		for (NatBitSet mec : mecs) {
			IntIterator iterator = mec.iterator();
			while (iterator.hasNext()) {
				int i = iterator.nextInt();
				// Check if any of the states in this MEC have not been touched yet - it is a new MEC then.
				if (!mecStates.contains(i)) {
					newMecs.add(mec);
					break;
				}
			}
		}

		if (verbose) {
			smgModelChecker.getLog().println(String.format("Found %d MECs, %d new ones", mecs.size(), newMecs.size()));
			smgModelChecker.getLog().flush();
		}

		for (NatBitSet mec : newMecs) {
			mecStates.or(mec);

			//Get all actions leaving MEC
			List<Object2DoubleMap<State>> actions = new ArrayList<>();
			mec.forEach((IntConsumer) state -> {
				for (int j = 0; j < partialMDP.getNumChoices(state); j++) {
					boolean successorsInSet = partialMDP.allSuccessorsInSet(state, j, mec);
					if (!successorsInSet) {
						actions.add(getDistribution(partialMDP, state, j));
					}
				}
			});

			if (!min && containsTarget(partialMDP, mec)) {
				// When the MEC contains the target, set all states of the MEC as target
				mec.forEach((IntConsumer) j -> {
					State state = partialMDP.getStatesList().get(j);
					stateUpdate.setQValue(state, new StateValue(1.0, 1.0));
					stateUpdate.setTarget(state, true);
				});
			} else {
				// When the MEC doesn't contain the target, the state's own actions are replaced by the outgoing transitions of the MEC
				// TODO: Better than this is to add a map to the cache which will point to the representative state
				if (!actions.isEmpty()) {
					mec.forEach((IntConsumer) j -> collapse(partialMDP, j, actions));
				}
			}
		}
		long duration = System.currentTimeMillis() - start;
		smgModelChecker.getLog().println(String.format("MEC collapsing done %s secs.", (double) duration / 1000));
		smgModelChecker.getLog().flush();
	}

	private boolean containsTarget(MDP mdp, NatBitSet mec) throws PrismException
	{
		IntIterator iterator = mec.iterator();
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			State state = mdp.getStatesList().get(i);
			if (stateUpdate.isTarget(state)) {
				return true;
			}
		}
		return false;
	}

	private void collapse(MDP mdp, int s, List<Object2DoubleMap<State>> actions)
	{
		CachedModelGenerator cmg = (CachedModelGenerator) modelGenerator;
		cmg.updateCache(mdp.getStatesList().get(s), actions);
	}

	private Object2DoubleMap<State> getDistribution(MDP mdp, int s, int i)
	{
		Object2DoubleMap<State> d = new Object2DoubleArrayMap<>();

		for (Int2DoubleMap.Entry e : mdp.getTransitions(s, i)) {
			d.put(mdp.getStatesList().get(e.getIntKey()), e.getDoubleValue());
		}

		return d;
	}
}