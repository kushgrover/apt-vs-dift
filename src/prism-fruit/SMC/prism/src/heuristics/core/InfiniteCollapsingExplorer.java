package heuristics.core;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import de.tum.in.naturals.unionfind.IntArrayUnionFind;
import de.tum.in.naturals.unionfind.IntUnionFind;
import explicit.Distribution;
import explicit.ECComputerFast;
import explicit.MDP;
import explicit.MDPExplicit;
import explicit.MDPSimple;
import explicit.MEC;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import prism.ModelGenerator;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.function.IntConsumer;
import java.util.stream.IntStream;

public class InfiniteCollapsingExplorer extends PrismComponent implements Explorer
{
	// TODO Sparse view on the partial model?
	// TODO Store the computed MECs directly?
	// TODO Predecessor-based optimizations (e.g. restricting the possible set of states investigated by MEC search or update of distributions after collapse)

	private final ModelExplorer partialExplorer;
	private final IntConsumer defaultAddedStateConsumer;
	private final MDPSimple collapsedModel = new MDPSimple();
	// UF structure keeping track of the representative state of all states in the partial model
	private final IntUnionFind mecUnionFind = new IntArrayUnionFind(0);
	private final NatBitSet statesInMEC = NatBitSets.set();
	private final NatBitSet collapsedExploredStates = NatBitSets.set();
	private boolean newStatesSinceLastCollapse = false;

	public InfiniteCollapsingExplorer(PrismComponent parent, ModelGenerator generator) throws PrismException
	{
		this(parent, generator, null);
	}

	public InfiniteCollapsingExplorer(PrismComponent parent, ModelGenerator generator, IntConsumer defaultAddedStateConsumer) throws PrismException
	{
		super(parent);
		this.defaultAddedStateConsumer = defaultAddedStateConsumer;
		partialExplorer = new ModelExplorer(this, generator, this::addState);

		IntIterator iterator = getPartialModel().getInitialStates().iterator();
		while (iterator.hasNext()) {
			int initialState = iterator.nextInt();
			exploreState(initialState);
			collapsedModel.addInitialState(initialState);
		}
	}

	public ModelExplorer getPartialExplorer()
	{
		return partialExplorer;
	}

	private void addState(int state)
	{
		int addedState = collapsedModel.addState();
		mecUnionFind.add();
		if (defaultAddedStateConsumer != null) {
			defaultAddedStateConsumer.accept(addedState);
		}
		assert addedState == state;
	}

	public Collection<NatBitSet> getPartialModelMECs()
	{
		Int2ObjectMap<NatBitSet> unionRootMap = new Int2ObjectAVLTreeMap<>();
		List<NatBitSet> foundBitSets = new ArrayList<>();
		IntIterator iterator = collapsedExploredStates.iterator();
		while (iterator.hasNext()) {
			int state = iterator.nextInt();
			int root = mecUnionFind.find(state);
			unionRootMap.computeIfAbsent(root, s -> {
				NatBitSet result = NatBitSets.set();
				foundBitSets.add(result);
				return result;
			}).set(state);
		}
		return foundBitSets;
	}

	@Override public void exploreState(int number) throws PrismException
	{
		assert collapsedModel.trans.size() > number : "State " + number + " was not added!";
		if (isStateExplored(number)) {
			return;
		}

		newStatesSinceLastCollapse = true;
		partialExplorer.exploreState(number);

		final MDP partialModel = partialExplorer.getPartialModel();
		final int numChoices = partialModel.getNumChoices(number);

		int addedActions = 0;
		for (int action = 0; action < numChoices; action++) {
			Distribution distribution = new Distribution();
			double outgoingWeight = 0d;
			// Collect all transitions which are not self loops. Since we are only concerned with infinite runs we can safely ignore them and reweigh all other
			// transitions
			for (Int2DoubleMap.Entry entry : partialModel.getTransitions(number, action)) {
				final int destination = mecUnionFind.find(entry.getIntKey());
				if (destination == number) {
					// Self-loop
					continue;
				}
				final double probability = entry.getDoubleValue();
				outgoingWeight += probability;
				distribution.add(destination, probability);
			}
			if (outgoingWeight == 0d) {
				// Only self-loops, completely ignore this action
				continue;
			}
			if (outgoingWeight < 1d) {
				// Some self-loops, rescale
				distribution.scale(outgoingWeight);
			}

			Object actionLabel = partialModel.getAction(number, action);
			if (actionLabel != null) {
				collapsedModel.addActionLabelledChoice(number, distribution, actionLabel);
			} else {
				collapsedModel.addChoice(number, distribution);
			}
			addedActions += 1;
		}

		if (addedActions == 0) {
			// The state had no actions or only self loops.
			collapsedModel.addDeadlockState(number);
			collapsedModel.trans.set(number, new Distribution[0]);
		}

		collapsedExploredStates.set(number);
		assert collapsedModel.trans.get(number) != null : numChoices;
	}

	@Override public boolean isStateExplored(int number)
	{
		return collapsedExploredStates.contains(number);
	}

	@Override public IntIterable getInitialStates()
	{
		return collapsedModel.getInitialStates();
	}

	@Override public int exploredStateCount()
	{
		return collapsedExploredStates.size();
	}

	public Int2ObjectMap<MEC> collapseMECs() throws PrismException
	{
		if (!newStatesSinceLastCollapse) {
			return new Int2ObjectArrayMap<>();
		}
		getLog().println("\nStarting MECs collapse", PrismLog.VL_HIGH);

		final List<MEC> mecs = ECComputerFast.computeMECs(collapsedModel);

		if (mecs.isEmpty()) {
			getLog().println("Found no MECs", PrismLog.VL_HIGH);
			return new Int2ObjectArrayMap<>();
		}
		NatBitSet statesInNewMECs = NatBitSets.setWithExpectedSize(collapsedModel.getNumStates());

		List<MEC> newMECs = new ArrayList<>();
		for (MEC mec : mecs) {
			NatBitSet states = mec.states;
			IntIterator iterator = states.iterator();
			while (iterator.hasNext()) {
				int state = iterator.nextInt();
				// Check if any of the states in this MEC have not been touched yet - it is a new MEC then.
				if (!statesInMEC.contains(state)) {
					statesInNewMECs.or(states);
					newMECs.add(mec);
					break;
				}
			}
		}
		statesInMEC.or(statesInNewMECs);
		if (getLog().isLoggable(PrismLog.VL_HIGH)) {
			getLog().println(String.format("Found %d new MECs with %d states", newMECs.size(), newMECs.stream().mapToInt(MEC::size).sum()), PrismLog.VL_HIGH);
		}

		// Collapse the states
		Int2ObjectMap<MEC> representativeMap = new Int2ObjectArrayMap<>();
		for (MEC mec : newMECs) {
			final NatBitSet states = mec.states;
			// Pick the representative for the new MEC
			int representativeState = mecUnionFind.find(states.firstInt());
			List<Distribution> distributions = new ArrayList<>(states.size() * 3);
			IntIterator iterator = states.iterator();
			while (iterator.hasNext()) {
				int state = iterator.nextInt();
				assert isStateExplored(state);
				// Union all states with the representative
				mecUnionFind.union(representativeState, state);
				// Delete the states transition from the collapsed model
				Distribution[] stateDistribution = collapsedModel.trans.set(state, null);
				distributions.addAll(Arrays.asList(stateDistribution));
			}
			// If the representative state we picked was part of a small zone, the above union may change the representative
			representativeState = mecUnionFind.find(representativeState);
			// States which are in a MEC but not the representative are counted as unexplored, the collapsed model does not contain them
			collapsedExploredStates.andNot(mec.states);
			collapsedExploredStates.set(representativeState);

			// Remove all transitions from the list which are completely inside of the MEC.
			distributions.removeIf(distribution -> distribution.isSubsetOf(states));
			collapsedModel.trans.set(representativeState, distributions.toArray(new Distribution[distributions.size()]));
			representativeMap.put(representativeState, mec);
		}

		// Check that the collapsedModel transition map is consistent with the exploredStates bit set.
		assert IntStream.range(0, collapsedModel.trans.size())
				.allMatch(state -> !collapsedExploredStates.contains(state) == (collapsedModel.trans.get(state) == null));

		if (getLog().isLoggable(PrismLog.VL_ALL)) {
			int transitionCount = 0;
			int actionCount = 0;
			int maxTransitions = 0;
			int maxActions = 0;
			int count = representativeMap.size();
			IntIterator iterator = representativeMap.keySet().iterator();
			while (iterator.hasNext()) {
				int representative = iterator.nextInt();
				Distribution[] distributions = collapsedModel.trans.get(representative);
				actionCount += distributions.length;
				maxActions = Math.max(maxActions, distributions.length);
				for (Distribution distribution : distributions) {
					transitionCount += distribution.size();
					maxTransitions = Math.max(maxTransitions, distribution.size());
				}
			}
			getLog().println(String.format("Collapsed states: %d, Actions: %.2f avg/%d max, Transitions %.2f avg/%d max",
					count, actionCount / (double) count, maxActions, transitionCount / (double) count, maxTransitions), PrismLog.VL_ALL);
		}

		// Remap all transitions. Other states might be pointing to some now merged state - we have to update them too.
		IntIterator iterator = collapsedExploredStates.iterator();
		while (iterator.hasNext()) {
			int state = iterator.nextInt();
			Distribution[] transitions = collapsedModel.trans.get(state);
			for (int action = 0; action < transitions.length; action++) {
				if (transitions[action].containsOneOf(statesInNewMECs)) {
					// Replace all successors by their corresponding UF-root.
					transitions[action] = transitions[action].map(mecUnionFind::find);
				} else {
					assert Objects.equals(transitions[action], transitions[action].map(mecUnionFind::find)) :
						String.format("%s vs %s", transitions[action], transitions[action].map(mecUnionFind::find));
				}
			}
		}

		return representativeMap;
	}

	public boolean isTransientState(int number)
	{
		return statesInMEC.contains(number);
	}

	@Override public Distribution getChoice(int state, int action)
	{
		assert isStateExplored(state);
		return collapsedModel.getChoice(state, action);
	}

	@Override public NatBitSet getExploredStates()
	{
		return collapsedExploredStates;
	}

	@Override public List<Distribution> getChoices(int state)
	{
		assert isStateExplored(state);
		Distribution[] choices = collapsedModel.getChoices(state);
		if (choices == null) {
			return Collections.emptyList();
		}
		return Arrays.asList(choices);
	}

	public MDPExplicit getPartialModel()
	{
		return partialExplorer.getPartialModel();
	}

	public MDPExplicit getCollapsedModel()
	{
		return collapsedModel;
	}

	public int partialModelExploredStateCount()
	{
		return partialExplorer.getExploredStates().size();
	}

	public int getCollapsedRepresentative(int stateNumber)
	{
		return mecUnionFind.find(stateNumber);
	}
}
