package heuristics.core;

import com.google.common.collect.Streams;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.Distribution;
import explicit.MDPExplicit;
import explicit.MEC;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntStack;
import prism.ModelGenerator;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismUtils;

import java.util.List;
import java.util.ListIterator;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.TimeUnit;

public class InfiniteCoreLearner extends PrismComponent
{
	private static final Random random = new Random();
	private static final int EXPLORES_PER_SAMPLE = 5;
	private static final int REPORT_PROGRESS_EVERY_STEPS = 50000;
	private static final boolean ACTION_PRUNING = true;

	private final InfiniteCollapsingExplorer explorer;
	private final IntList visitedStates = new IntArrayList();
	// Can also be implemented as array, but this enables a sparse representation, saving memory
	private final Int2DoubleMap reachabilityUpperBounds = new Int2DoubleLinkedOpenHashMap();
	private final Int2ObjectMap<NatBitSet> solvedActionsMap;
	private final boolean enableActionPruning;
	private double precision = 1;
	private int sampleSteps = 0;
	private int samples = 0;
	private int backtraceCount = 0;
	private int backtraceSteps = 0;
	private int loopCount = 0;
	private int collapseThreshold = 1;

	public InfiniteCoreLearner(PrismComponent parent, ModelGenerator model, boolean enableActionPruning) throws PrismException
	{
		super(parent);
		this.explorer = new InfiniteCollapsingExplorer(this, model, null);
		reachabilityUpperBounds.defaultReturnValue(1);
		this.enableActionPruning = enableActionPruning;
		solvedActionsMap = enableActionPruning ? new Int2ObjectLinkedOpenHashMap<>() : null;
	}

	public InfiniteCoreLearner(PrismComponent parent, ModelGenerator model) throws PrismException
	{
		this(parent, model, ACTION_PRUNING);
	}

	public InfiniteCollapsingExplorer getCollapsingExplorer()
	{
		return explorer;
	}

	private boolean same(double first, double second)
	{
		assert precision >= 0d;
		return Math.abs(first - second) < precision;
	}

	private int getBestAction(int state)
	{
		// Note: We don't call getBestAction on states with reachability upper bound less than precision
		assert explorer.isStateExplored(state) : String.format("State %d not explored! Current path: %s", state, visitedStates);
		assert reachabilityUpperBounds.get(state) >= precision;

		List<Distribution> choices = explorer.getChoices(state);
		if (choices.isEmpty()) {
			return -1;
		}
		if (choices.size() == 1) {
			return 0;
		}

		double bestValue = 0d;
		NatBitSet solvedActions = enableActionPruning ? solvedActionsMap.get(state) : null;
		// Over-approximate the set of precision-best actions
		final IntList bestActions = new IntArrayList();
		final DoubleList bestActionValues = new DoubleArrayList(choices.size());
		for (ListIterator<Distribution> iterator = choices.listIterator(); iterator.hasNext(); ) {
			// Iterate over all successor distributions and gather their values
			final int action = iterator.nextIndex();
			final Distribution distribution = iterator.next();
			if (solvedActions != null && solvedActions.contains(action)) {
				continue;
			}
			final double value = distribution.sumWeighted(reachabilityUpperBounds);
			if (value > bestValue) {
				// This value will be the new best value
				if (value - bestValue > precision) {
					// This value actually is significantly better than any of the stored values
					bestActions.clear();
					bestActionValues.clear();
				}
				bestActions.add(action);
				bestActionValues.add(value);
				bestValue = value;
			} else if (value > bestValue - precision) {
				// Its not a new best but good enough (within (bestVal - precision, bestVal) interval)
				bestActions.add(action);
				bestActionValues.add(value);
			}
		}
		// This would fail if action pruning is enabled and the UB is less than precision, since then all successor actions necessarily are marked as solved and
		// no action is added.
		assert !bestActions.isEmpty() : "No best action available";
		if (bestActions.size() == 1) {
			// No need to randomize
			return bestActions.getInt(0);
		}

		// Properly filter the set now - the bestAction list might contain actions which are n*precision worse than the best value.
		// This list will contain the indices of the bestActions list.
		final IntList bestActionIndices = new IntArrayList(bestActionValues.size());
		for (int i = 0, size = bestActionValues.size(); i < size; i++) {
			double value = bestActionValues.getDouble(i);
			assert value <= bestValue;
			if (same(value, bestValue)) {
				bestActionIndices.add(i);
			}
		}
		// There has to be a witness for the bestValue
		assert !bestActionIndices.isEmpty();

		// If we only have one best value there is no need to randomize
		final int bestActionIndex = bestActionIndices.size() == 1 ? 0 : random.nextInt(bestActionIndices.size());
		return bestActions.getInt(bestActionIndices.getInt(bestActionIndex));
	}

	private int sampleNextState(int state, int action)
	{
		assert explorer.isStateExplored(state);
		final Distribution distribution = explorer.getChoice(state, action);

		// Sample according to max reachability upper bound
		final int[] successors = distribution.getSupport();
		assert successors.length > 0;
		if (successors.length == 1) {
			final int potentialSuccessor = successors[0];
			if (reachabilityUpperBounds.get(potentialSuccessor) < precision) {
				// Always return -1 if there is nothing worth exploring - this successor is already solved. Its values will be propagated onto the current
				// state by update
				return -1;
			}
			return potentialSuccessor;
		}

		// Get sum of all successors' upper bounds
		double sum = 0d;
		double[] values = new double[successors.length];
		for (int i = 0; i < successors.length; i++) {
			final double value = reachabilityUpperBounds.get(successors[i]);
			if (value < precision) {
				// Don't consider already solved states
				continue;
			}
			sum += value;
			values[i] = value;
		}
		if (sum == 0d) {
			// Only solved successors, again return -1 to signal this
			return -1;
		}
		assert precision <= sum && sum <= successors.length;

		// Sample a random value in [0, sum)
		final double sampledValue = random.nextDouble() * sum;
		// Search the successor corresponding to this value
		double partialSum = 0d;
		for (int i = 0; i < successors.length; i++) {
			final double value = values[i];
			if (value == 0d) {
				continue;
			}
			assert value >= precision;
			partialSum += value;
			if (partialSum >= sampledValue) {
				return successors[i];
			}
		}
		throw new IllegalStateException("Not sampling any state");
	}

	public MDPExplicit learn(double precision) throws PrismException
	{
		if (this.precision < precision) {
			assert Streams.stream(explorer.getInitialStates())
					.mapToInt(explorer::getCollapsedRepresentative)
					.distinct()
					.allMatch(representative -> reachabilityUpperBounds.get(representative) < precision);
			if (getLog().isLoggable(PrismLog.VL_HIGH)) {
				getLog().println(String.format("Returning already computed model of precision %f (requested: %f)", this.precision, precision),
						PrismLog.VL_HIGH);
			}
			return explorer.getPartialModel();
		}
		getLog().println(String.format("Learning the infinite core for precision %g from initial states %s", precision, explorer.getInitialStates()));
		this.precision = precision;

		long timer = System.nanoTime();

		IntIterator iterator = explorer.getInitialStates().iterator();
		while (iterator.hasNext()) {
			int initialState = iterator.nextInt();
			// The representative of the initial states might be a different state
			int representative = explorer.getCollapsedRepresentative(initialState);
			while (reachabilityUpperBounds.get(representative) > this.precision) {
				if (exploreAndUpdate(representative)) {
					// If MECs have been merged, update the representative (it might have changed)
					representative = explorer.getCollapsedRepresentative(initialState);
					assert explorer.getCollapsedRepresentative(representative) == explorer.getCollapsedRepresentative(initialState)
							&& explorer.isStateExplored(representative);
				}
			}
		}
		explorer.collapseMECs();

		timer = System.nanoTime() - timer;

		if (getLog().isLoggable(PrismLog.VL_DEFAULT)) {
			final String progressString = String.format("%n== Finished core learning (precision %g) ==%n"
							+ "  Total samples: %d, steps: %d, avg len: %f, backtrace count: %d, steps: %d%n"
							+ "  States: %d/%d in partial/collapsed model%n"
							+ "  Time: %f sec%n",
					precision, samples, sampleSteps, (double) sampleSteps / (double) samples, backtraceCount, backtraceSteps,
					explorer.partialModelExploredStateCount(), explorer.exploredStateCount(),
					timer / (double) TimeUnit.SECONDS.toNanos(1));
			getLog().println(progressString, PrismLog.VL_DEFAULT);
		}

		return explorer.getPartialModel();
	}

	private boolean exploreAndUpdate(int state) throws PrismException
	{
		assert visitedStates.isEmpty();
		assert reachabilityUpperBounds.get(state) >= precision;
		samples++;

		IntStack visitStack = (IntStack) visitedStates;
		int currentState = state;
		int exploreCount = 0;

		// Sample a path
		while (!visitedStates.contains(currentState)) {
			sampleSteps++;
			if (sampleSteps % REPORT_PROGRESS_EVERY_STEPS == 0) {
				reportProgress();
			}

			// Determine the next best action
			int bestAction = getBestAction(currentState);
			if (bestAction == -1) {
				// No successor available
				reachabilityUpperBounds.put(currentState, 0d);
				break;
			}
			visitStack.push(currentState);
			// Sample the successor
			currentState = sampleNextState(currentState, bestAction);
			assert currentState == -1 || reachabilityUpperBounds.get(currentState) >= precision;
			if (currentState == -1) {
				// We won't find anything of value if we continue to follow this path, backtrace until we find an interesting state again
				// Note: We might as well completely restart the sampling here, but then we potentially have to move to this interesting "fringe" region again
				do {
					backtraceSteps++;
					currentState = visitStack.popInt();
				} while (update(currentState) < precision);
				backtraceCount++;
			} else if (!explorer.isStateExplored(currentState)) {
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

		// Collapse MECs
		if (visitedStates.contains(currentState)) {
			// Sampled the same state twice
			assert explorer.isStateExplored(currentState);
			loopCount++;
			// We looped quite often - chances for this are high if there is a MEC, otherwise the sampling probabilities would decrease
			if (loopCount > collapseThreshold) {
				// Search for MECs in the partial MDP and process the MECs which were not processed yet
				updateMecs(explorer.collapseMECs());

				loopCount = 0;
				collapseThreshold = explorer.exploredStateCount();

				// The path could now contain invalid states (i.e. some which have been merged), we have to sample again.
				visitedStates.clear();
				return true;
			}
		}

		// Propagate values backwards along the path
		while (!visitStack.isEmpty()) {
			int visitedState = visitStack.popInt();
			update(visitedState);
		}
		return false;
	}

	public void reportProgress()
	{
		if (!getLog().isLoggable(PrismLog.VL_HIGH)) {
			return;
		}
		Int2DoubleMap initialStateValues = new Int2DoubleArrayMap();
		explorer.getPartialModel().getInitialStates().forEach((int initialState) -> {
			int collapsedRepresentative = explorer.getCollapsedRepresentative(initialState);
			initialStateValues.put(initialState, reachabilityUpperBounds.get(collapsedRepresentative));
		});

		String progressString = String.format("%n== Progress Report ==%n"
						+ "  Trials: %d, steps: %d, avg len: %f, backtrace count: %d, steps: %d%n"
						+ "  States: %d/%d in partial/collapsed model%n"
						+ "  Bounds in initial states: %s",
				samples, sampleSteps, (double) sampleSteps / (double) samples, backtraceCount, backtraceSteps,
				explorer.partialModelExploredStateCount(), explorer.exploredStateCount(),
				initialStateValues);
		if (enableActionPruning) {
			progressString += String.format("%n  Pruned actions: %d in %d states",
					solvedActionsMap.values().stream().mapToInt(Set::size).sum(), solvedActionsMap.size());
		}
		getLog().println(progressString, PrismLog.VL_HIGH);
	}

	private void updateMecs(Int2ObjectMap<MEC> mecs)
	{
		if (mecs.isEmpty()) {
			return;
		}

		for (Int2ObjectMap.Entry<MEC> entry : mecs.int2ObjectEntrySet()) {
			final int representative = entry.getIntKey();
			assert explorer.isStateExplored(representative);
			if (enableActionPruning) {
				// The representatives actions have changed, drop the cache
				solvedActionsMap.remove(representative);
			}
			// Update the representative - it represents the whole MEC now!
			update(representative);

			// This is only to save memory
			final NatBitSet states = entry.getValue().states;
			IntIterator iterator = states.iterator();
			while (iterator.hasNext()) {
				int i = iterator.nextInt();
				if (i == representative) {
					continue;
				}
				if (enableActionPruning) {
					solvedActionsMap.remove(i);
				}
				reachabilityUpperBounds.remove(i);
			}
		}
	}

	private double update(int state)
	{
		// Determine the maximal expected upper bound of the successors
		double maximalValue = 0d;

		boolean newSetCreated = false;
		NatBitSet solvedActions = enableActionPruning ? solvedActionsMap.get(state) : null;
		final List<Distribution> choices = explorer.getChoices(state);
		for (ListIterator<Distribution> iterator = choices.listIterator(); iterator.hasNext(); ) {
			int action = iterator.nextIndex();
			Distribution distribution = iterator.next();
			if (solvedActions != null && solvedActions.contains(action)) {
				assert distribution.sumWeighted(reachabilityUpperBounds) < precision;
				continue;
			}

			final double actionValue = distribution.sumWeighted(reachabilityUpperBounds);
			assert 0 <= actionValue && (actionValue <= 1d || PrismUtils.doublesAreEqual(actionValue, 1d)) : String.format("Invalid value %.20f", actionValue);
			if (enableActionPruning && actionValue < precision) {
				if (solvedActions == null) {
					solvedActions = NatBitSets.boundedSet(choices.size());
					newSetCreated = true;
				}
				solvedActions.set(action);
			}
			maximalValue = Math.max(actionValue, maximalValue);
		}

		assert 0d <= maximalValue && (maximalValue <= 1d || PrismUtils.doublesAreEqual(maximalValue, 1d)) : String.format("Invalid value %.20f", maximalValue);

		if (enableActionPruning) {
			if (maximalValue < precision) {
				// We won't visit this state anymore
				if (!newSetCreated) {
					solvedActionsMap.remove(state);
				}
			} else if (newSetCreated) {
				solvedActionsMap.put(state, solvedActions);
			}
		}

		if (maximalValue >= 1d) {
			assert reachabilityUpperBounds.get(state) == 1d : String.format("Maximal value is 1 but bound was %f", reachabilityUpperBounds.get(state));
			// 1 is the default value, no need to store it.
			// reachabilityUpperBounds.remove(state);
			return 1d;
		}
		reachabilityUpperBounds.put(state, maximalValue);
		return maximalValue;
	}
}
