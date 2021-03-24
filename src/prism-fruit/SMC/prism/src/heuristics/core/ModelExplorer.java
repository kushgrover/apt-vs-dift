package heuristics.core;

import com.google.common.collect.Iterables;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.Distribution;
import explicit.MDPExplicit;
import explicit.MDPSimple;
import it.unimi.dsi.fastutil.ints.Int2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.objects.Object2IntAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import parser.State;
import prism.ModelGenerator;
import prism.PrismComponent;
import prism.PrismException;

import java.util.Arrays;
import java.util.List;
import java.util.function.IntConsumer;

public class ModelExplorer extends PrismComponent implements Explorer
{
	private final MDPSimple partialModel = new MDPSimple();
	private final ModelGenerator generator;
	private final IntConsumer defaultAddedStateConsumer;
	private final Object2IntMap<State> stateMap = new Object2IntAVLTreeMap<>();
	private final Int2ObjectMap<State> indexMap = new Int2ObjectAVLTreeMap<>();
	// All states which are in the partial model and explored
	private final NatBitSet exploredStates = NatBitSets.set();

	public ModelExplorer(PrismComponent parent, ModelGenerator generator, IntConsumer defaultAddedStateConsumer) throws PrismException
	{
		super(parent);
		this.generator = generator;
		this.defaultAddedStateConsumer = defaultAddedStateConsumer;
		stateMap.defaultReturnValue(-1);

		for (State initialState : generator.getInitialStates()) {
			int stateNumber = addState(initialState, defaultAddedStateConsumer);
			exploreState(stateNumber);
			partialModel.addInitialState(stateNumber);
		}
	}

	public ModelGenerator getGenerator()
	{
		return generator;
	}

	@Override public NatBitSet getExploredStates()
	{
		return exploredStates;
	}

	protected MDPExplicit getPartialModel()
	{
		return partialModel;
	}

	void exploreState(int stateNumber, IntConsumer addedStateConsumer) throws PrismException
	{
		assert indexMap.containsKey(stateNumber) && stateMap.containsKey(indexMap.get(stateNumber));
		if (isStateExplored(stateNumber)) {
			return;
		}
		exploredStates.set(stateNumber);

		State state = indexMap.get(stateNumber);
		assert state != null;
		generator.exploreState(state);

		int actions = generator.getNumChoices();
		int totalTransitions = 0;
		for (int action = 0; action < actions; action++) {
			int transitionCount = generator.getNumTransitions(action);
			totalTransitions += transitionCount;

			Distribution distribution = new Distribution();
			for (int transition = 0; transition < transitionCount; transition++) {
				final double transitionProbability = generator.getTransitionProbability(action, transition);
				final State transitionTarget = generator.computeTransitionTarget(action, transition);
				final int transitionTargetNumber = addState(transitionTarget, addedStateConsumer);
				distribution.add(transitionTargetNumber, transitionProbability);
			}
			final Object actionLabel = generator.getChoiceAction(action);
			final double stateActionReward = generator.getStateActionReward(0, state, actionLabel);
			int retVal = partialModel.addActionLabelledChoice(stateNumber, distribution, actionLabel);
			assert retVal != -1;
		}
		assert partialModel.getNumChoices(stateNumber) == actions;
		if (totalTransitions == 0) {
			partialModel.addDeadlockState(stateNumber);
		}
	}

	@Override public void exploreState(int stateNumber) throws PrismException
	{
		exploreState(stateNumber, defaultAddedStateConsumer);
	}

	@Override public boolean isStateExplored(int number)
	{
		return exploredStates.contains(number);
	}

	private int addState(State state, IntConsumer addedStateConsumer)
	{
		assert state != null;
		int stateNumber = stateMap.getInt(state);
		if (stateNumber != -1) {
			return stateNumber;
		}
		int newStateNumber = partialModel.addState();
		assert newStateNumber == stateMap.size();
		stateMap.put(state, newStateNumber);
		indexMap.put(newStateNumber, state);
		if (addedStateConsumer != null) {
			addedStateConsumer.accept(newStateNumber);
		}
		return newStateNumber;
	}

	@Override public Distribution getChoice(int state, int action)
	{
		assert isStateExplored(state);
		return partialModel.getChoice(state, action);
	}

	@Override public int getInitialState()
	{
		return Iterables.getOnlyElement(partialModel.getInitialStates());
	}

	@Override public IntIterable getInitialStates()
	{
		return partialModel.getInitialStates();
	}

	@Override public int exploredStateCount()
	{
		return exploredStates.size();
	}

	@Override public List<Distribution> getChoices(int state)
	{
		assert isStateExplored(state);
		return Arrays.asList(partialModel.getChoices(state));
	}

	public State getState(int stateNumber)
	{
		return indexMap.get(stateNumber);
	}
}
