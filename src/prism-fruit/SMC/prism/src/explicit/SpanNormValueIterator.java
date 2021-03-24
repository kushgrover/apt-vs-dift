package explicit;

import common.Time;
import it.unimi.dsi.fastutil.ints.Int2DoubleFunction;
import it.unimi.dsi.fastutil.ints.Int2DoubleLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntHash;
import it.unimi.dsi.fastutil.ints.IntIterator;
import prism.PrismComponent;
import prism.PrismSettings;
import prism.PrismUtils;

import java.util.function.IntPredicate;
import java.util.function.IntToDoubleFunction;

public abstract class SpanNormValueIterator extends PrismComponent
{
	private static class IdentityHash implements IntHash.Strategy {

		@Override public int hashCode(int e)
		{
			return e;
		}

		@Override public boolean equals(int a, int b)
		{
			return a == b;
		}
	}

	final boolean useAdaptiveTau;
	final int adaptiveTauTransformationCheckSteps;
	final double adaptiveTauTargetValue;
	Int2DoubleFunction currentValues;
	Int2DoubleFunction previousValues;

	SpanNormValueIterator(PrismComponent parent)
	{
		super(parent);
		this.useAdaptiveTau = settings.getBoolean(PrismSettings.PRISM_MDP_MP_ADAPTIVE_TAU);
		this.adaptiveTauTransformationCheckSteps = settings.getInteger(PrismSettings.PRISM_MDP_MP_ADAPTIVE_TAU_STEPS);
		this.adaptiveTauTargetValue = settings.getDouble(PrismSettings.PRISM_MDP_MP_APERIODICITYTAU);
		// TODO Improve this - factor 2 faster!
		// Solution: Delegate storage accesses to the runner, pull run to the Iterator and expose all the variables?
		currentValues = new Int2DoubleLinkedOpenHashMap();
		previousValues = new Int2DoubleLinkedOpenHashMap();
	}

	abstract Model getModel();

	private void swap()
	{
		Int2DoubleFunction swap = currentValues;
		currentValues = previousValues;
		previousValues = swap;
	}

	abstract class Runner
	{
		// Model / problem fields
		private final IntPredicate checkStepPredicate;
		private final double precision;

		// Progress fields
		private int iterations = 0;
		private boolean adaptiveTau;
		private double currentTauValue;
		private double lastAdaptiveTauCheckSpanNorm;
		private int lastAdaptiveTauCheck;
		private double previousSpanNorm;
		private double currentSpanNorm;
		private boolean converged;

		// Result fields
		private double resultMaximumValue = Double.NaN;
		private double resultMinimumValue = Double.NaN;
		private double timeTaken = 0d;
		private double resultMeanPayoff = Double.NaN;

		Runner(IntPredicate checkStepPredicate, double precision)
		{
			assert precision >= 0.0d;
			this.checkStepPredicate = checkStepPredicate;
			this.precision = precision;

			// Do we choose the tau adaptively
			adaptiveTau = useAdaptiveTau;
			// If tau is adaptive, we start with 1 and switch if needed, otherwise we start with the stored tau value
			currentTauValue = useAdaptiveTau ? 1d : adaptiveTauTargetValue;
			// This stores the span norm of the last time we checked if we need to enable the transformation
			lastAdaptiveTauCheckSpanNorm = Double.POSITIVE_INFINITY;
			lastAdaptiveTauCheck = 0;

			// Keep track of the span norm
			previousSpanNorm = Double.POSITIVE_INFINITY;
			currentSpanNorm = Double.NaN;
			converged = false;
		}

		abstract double operator(int state, IntToDoubleFunction stateValues, double tau);

		public double getPrecision()
		{
			return precision;
		}

		public int getIterations()
		{
			return iterations;
		}

		public boolean isConverged()
		{
			return converged;
		}

		public int getFirstState() {
			return getStateIterator().nextInt();
		}

		public double getResultMaximumValue()
		{
			assert !Double.isNaN(resultMaximumValue) && isConverged();
			return resultMaximumValue;
		}

		public double getResultMinimumValue()
		{
			assert !Double.isNaN(resultMinimumValue) && isConverged();
			return resultMinimumValue;
		}

		public double getTimeTaken()
		{
			return timeTaken;
		}

		abstract IntIterator getStateIterator();

		public double getResultMeanPayoff()
		{
			assert !Double.isNaN(resultMeanPayoff) && isConverged();
			return resultMeanPayoff;
		}

		public double getValue(int state)
		{
			return currentValues.get(state) / currentTauValue;
		}

		public double getPreviousValue(int state)
		{
			return previousValues.get(state) / currentTauValue;
		}

		public double getCurrentGain(int state) {
			return (currentValues.get(state) - previousValues.get(state)) / currentTauValue;
		}

		public double getObtainedPrecision()
		{
			assert isConverged();
			return currentSpanNorm;
		}

		public boolean run(int steps)
		{
			if (converged) {
				return true;
			}

			// How many steps to do at most
			long stepBound = steps >= 0 ? steps : Long.MAX_VALUE;
			Time time = new Time();
			int firstState = getFirstState();

			for (long step = 0; step < stepBound; step++) {
				if (step > 0 && checkStepPredicate.test(iterations)) {
					// This is a step where we should check for convergence

					double minimumDifference = Double.POSITIVE_INFINITY;
					double maximumDifference = Double.NEGATIVE_INFINITY;

					swap();
					IntIterator iterator = getStateIterator();
					while (iterator.hasNext()) {
						// Iterate over all states and apply the operator on each state
						int state = iterator.nextInt();
						double value = operator(state, previousValues::get, currentTauValue);
						currentValues.put(state, value);

						double difference = value - previousValues.get(state);
						assert difference >= -1e-8 : String.format("Different values or overflow: current %s, previous %s in iteration %d",
								currentValues.get(state), previousValues.get(state), iterations);
						if (difference > maximumDifference) {
							maximumDifference = difference;
						}
						if (difference < minimumDifference) {
							minimumDifference = difference;
						}
					}
					assert Double.isFinite(minimumDifference) && Double.isFinite(maximumDifference) && minimumDifference <= maximumDifference;
					currentSpanNorm = maximumDifference - minimumDifference;
					assert currentSpanNorm >= 0d;
					assert previousSpanNorm >= currentSpanNorm
							|| PrismUtils.doublesAreCloseAbs(previousSpanNorm, currentSpanNorm, precision * 2d) :
							String.format("Different values or overflow: %g -> %g", previousSpanNorm, currentSpanNorm);

					// If span norm of the difference of the current and previous iteration is less than epsilon, stop
					if (currentSpanNorm < currentTauValue * precision) {
						converged = true;
						break;
					}
					previousSpanNorm = currentSpanNorm;

					if (adaptiveTau && iterations - lastAdaptiveTauCheck > adaptiveTauTransformationCheckSteps) {
						if (lastAdaptiveTauCheckSpanNorm - currentSpanNorm < precision) {
							// Span norm didn't change much - decrease tau value
							currentTauValue = adaptiveTauTargetValue;
							adaptiveTau = false;

							// Adapt the state values
							IntIterator checkStates = getStateIterator();
							while (checkStates.hasNext()) {
								int state = checkStates.nextInt();
								currentValues.put(state, currentValues.get(state) * currentTauValue);
								previousValues.put(state, previousValues.get(state) * currentTauValue);
							}
							// Reset the span norm checking
							previousSpanNorm = Double.POSITIVE_INFINITY;
						}
						lastAdaptiveTauCheckSpanNorm = currentSpanNorm;
						lastAdaptiveTauCheck = iterations;
					}
				} else {
					// Just perform VI without checking the values
					swap();
					IntIterator iterator = getStateIterator();
					while (iterator.hasNext()) {
						int state = iterator.nextInt();
						double value = operator(state, previousValues::get, currentTauValue);
						currentValues.put(state, value);
					}
				}

				// Relative VI: Re-normalize the values to avoid numerical errors on large models
				if (currentValues.get(firstState) > 1000d) {
					// TODO Make this controllable, the 1000d switch is quite arbitrary
					double referenceValue = previousValues.get(firstState);
					IntIterator iterator = getStateIterator();
					while (iterator.hasNext()) {
						int state = iterator.nextInt();
						currentValues.put(state, currentValues.get(state) - referenceValue);
						previousValues.put(state, previousValues.get(state) - referenceValue);
					}
				}
				iterations += 1;
			}

			if (!converged) {
				return false;
			}

			// Compute relative differences to avoid numerical errors and overflows
			double relativeValue = (currentValues.get(firstState) - previousValues.get(firstState)) / currentTauValue;
			double accumulatedValue = 0d;
			double maximumValue = relativeValue;
			double minimumValue = relativeValue;

			int stateCount = 0;
			IntIterator iterator = getStateIterator();
			while (iterator.hasNext()) {
				int state = iterator.nextInt();
				double currentValue = currentValues.get(state);
				double difference = (currentValue - previousValues.get(state)) / currentTauValue;
				assert difference >= -PrismUtils.epsilonDouble;
				maximumValue = Math.max(maximumValue, difference);
				minimumValue = Math.min(minimumValue, difference);
				accumulatedValue += difference - relativeValue;
				stateCount += 1;
			}

			// Compute the average reward by re-weighting with tau and using the fact that v_n approximates n*gain + bias
			resultMeanPayoff = (relativeValue + (accumulatedValue / (double) stateCount));
			assert !Double.isNaN(resultMeanPayoff);

			resultMaximumValue = maximumValue;
			resultMinimumValue = minimumValue;
			timeTaken += time.elapsedSeconds();
			return true;
		}

	}
}
