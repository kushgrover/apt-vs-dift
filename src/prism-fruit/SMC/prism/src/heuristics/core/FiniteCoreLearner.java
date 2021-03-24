package heuristics.core;

import explicit.Distribution;
import explicit.MDPExplicit;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.Int2DoubleAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntStack;
import prism.ModelGenerator;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.function.IntConsumer;

public class FiniteCoreLearner extends PrismComponent
{
	private static final Random random = new Random();
	private static final int EXPLORES_PER_SAMPLE = 5;
	private static final int REPORT_PROGRESS_EVERY_STEPS = 50000;

	private final ModelExplorer explorer;
	private final IntList visitedStates = new IntArrayList();
	private final List<Int2DoubleMap> reachabilityUpperBounds = new ArrayList<>();
	private int sampleSteps = 0;
	private int samples = 0;
	private double precision = 1;
	private int stepBound;

	public FiniteCoreLearner(PrismComponent parent, ModelGenerator generator) throws PrismException
	{
		super(parent);
		this.explorer = new ModelExplorer(this, generator, this::addState);
	}

	private void addState(int stateNumber)
	{
		assert reachabilityUpperBounds.size() == stateNumber;
		Int2DoubleMap map = new Int2DoubleAVLTreeMap();
		map.defaultReturnValue(1);
		reachabilityUpperBounds.add(map);
	}

	public MDPExplicit learn(int stepBound, double precision) throws PrismException
	{
		getLog().println(String.format("Learning the infinite core for step bound %d and precision %g", stepBound, precision));
		this.precision = precision;
		this.stepBound = stepBound;

		long timer = System.nanoTime();

		IntIterator iterator = explorer.getInitialStates().iterator();
		while (iterator.hasNext()) {
			int initialState = iterator.nextInt();
			while (getValue(initialState, 0) > this.precision) {
				exploreAndUpdate(initialState);
			}
		}

		timer = System.nanoTime() - timer;

		if (getLog().isLoggable(PrismLog.VL_HIGH)) {
			final String progressString = String.format("%n== Finished core learning (precision %g) ==%n"
							+ "  Total samples: %d, steps: %d, avg len: %f%n"
							+ "  States: %d in partial model%n"
							+ "  Time: %f sec%n",
					precision, samples, sampleSteps, (double) sampleSteps / (double) samples,
					explorer.exploredStateCount(),
					timer / (double) TimeUnit.SECONDS.toNanos(1));
			getLog().println(progressString, PrismLog.VL_HIGH);
		}

		return explorer.getPartialModel();
	}

	private boolean same(double first, double second)
	{
		assert precision >= 0d;
		return Math.abs(first - second) < precision;
	}

	private int getBestAction(int state, int remainingSteps)
	{
		assert explorer.isStateExplored(state) && remainingSteps > 0;
		List<Distribution> choices = explorer.getChoices(state);
		if (choices.isEmpty()) {
			return -1;
		}
		if (choices.size() == 1) {
			return 0;
		}

		double bestValue = 0d;
		DoubleList actionValues = new DoubleArrayList(choices.size());
		for (Distribution distribution : choices) {
			double value = distribution.sumWeighted(s -> getValue(s, remainingSteps - 1));
			actionValues.add(value);
			if (value > bestValue) {
				bestValue = value;
			}
		}
		IntList bestActions = new IntArrayList();
		for (int i = 0; i < actionValues.size(); i++) {
			double value = actionValues.getDouble(i);
			assert value <= bestValue;
			if (same(value, bestValue)) {
				bestActions.add(i);
			}
		}
		assert !bestActions.isEmpty();
		if (bestActions.size() == 1) {
			return bestActions.getInt(0);
		}
		return bestActions.getInt(random.nextInt(bestActions.size()));
	}

	private int sampleNextState(int state, int action, int remainingSteps)
	{
		assert explorer.isStateExplored(state) && remainingSteps > 0;
		final Distribution distribution = explorer.getChoice(state, action);

		// Sample according to max reachability upper bound
		final int[] successors = distribution.getSupport();
		assert successors.length > 0;
		if (successors.length == 1) {
			return successors[0];
		}

		// Get sum of all successors
		double sum = 0d;
		double[] values = new double[successors.length];
		for (int i = 0; i < successors.length; i++) {
			final double value = getValue(successors[i], remainingSteps - 1);
			if (value < precision) {
				// Don't move to already solved states
				continue;
			}
			sum += value;
			values[i] = value;
		}

		assert precision < sum && sum <= successors.length;

		// Sample a random value in [0, sum)
		final double sampledValue = random.nextDouble() * sum;
		// Search the successor corresponding to this value
		double partialSum = 0d;
		for (int i = 0; i < successors.length; i++) {
			partialSum += values[i];
			if (partialSum >= sampledValue) {
				return successors[i];
			}
		}
		throw new IllegalStateException("Not sampling any state");
	}

	private void exploreAndUpdate(int state) throws PrismException
	{
		assert visitedStates.isEmpty();
		samples++;

		IntStack visitStack = (IntStack) visitedStates;
		int currentState = state;
		int remainingSteps = stepBound;
		int exploreCount = 0;

		// Sample a path
		while (true) {
			sampleSteps++;
			if (sampleSteps % REPORT_PROGRESS_EVERY_STEPS == 0) {
				reportProgress();
			}

			// Determine the next best action
			int bestAction = getBestAction(currentState, remainingSteps);
			if (bestAction == -1) {
				// No successor available
				setValue(currentState, remainingSteps, 0d);
				break;
			}
			visitStack.push(currentState);
			if (remainingSteps == 1) {
				// No need to sample a successor - remaining steps = 0 means value = 0 anyway
				break;
			}
			// Sample the successor
			currentState = sampleNextState(currentState, bestAction, remainingSteps);
			remainingSteps -= 1;
			// We won't find anything of value if we continue to follow this path
			if (currentState == -1) {
				break;
			}

			if (!explorer.isStateExplored(currentState)) {
				// Found a new state

				if (exploreCount == EXPLORES_PER_SAMPLE) {
					// Explored along this path long enough
					break;
				}

				// Explore the state and continue sampling from here
				exploreCount++;
				explorer.exploreState(currentState);
			}
		}

		// Propagate values backwards along the path
		while (!visitStack.isEmpty()) {
			assert stepBound - visitedStates.size() == remainingSteps;
			int visitedState = visitStack.popInt();
			update(visitedState, remainingSteps);
			remainingSteps += 1;
		}
	}

	public void reportProgress()
	{
		if (!getLog().isLoggable(PrismLog.VL_HIGH)) {
			return;
		}
		Int2DoubleMap initialStateValues = new Int2DoubleArrayMap();
		explorer.getPartialModel().getInitialStates().forEach((IntConsumer) initialState -> initialStateValues.put(initialState, getValue(initialState, 0)));

		final String progressString = String.format("%n== Progress Report ==%n"
						+ "  Trials: %d, steps: %d, avg len: %f%n"
						+ "  States: %d in partia model%n"
						+ "  Bounds in initial states: %s",
				samples, sampleSteps, (double) sampleSteps / (double) samples,
				explorer.exploredStateCount(),
				initialStateValues);
		getLog().println(progressString, PrismLog.VL_HIGH);
	}

	private void setValue(int state, int remainingSteps, double value) {
		assert 0 < remainingSteps && 0 <= state && state < reachabilityUpperBounds.size();
		final Int2DoubleMap map = reachabilityUpperBounds.get(state);
		if (value == 1d) {
			map.remove(remainingSteps);
		} else {
			map.put(remainingSteps, value);
		}
	}

	private double getValue(int state, int remainingSteps)
	{
		assert 0 <= remainingSteps && 0 <= state && state < reachabilityUpperBounds.size();
		if (remainingSteps == 0) {
			return 0d;
		}
		return reachabilityUpperBounds.get(state).get(remainingSteps);
	}

	private void update(int state, int remainingSteps)
	{
		// Determine the maximal expected upper bound of the successors
		double maximalValue = 0d;
		for (Distribution distribution : explorer.getChoices(state)) {
			maximalValue = Math.max(distribution.sumWeighted(s -> getValue(s, remainingSteps - 1)), maximalValue);
		}
		assert 0d <= maximalValue && maximalValue <= 1d;
		setValue(state, remainingSteps, maximalValue);
	}
}

