package explicit.result;

import java.util.function.IntToDoubleFunction;

public class IntervalResultSparse implements IntervalResult
{
	private final IntToDoubleFunction upper;
	private final IntToDoubleFunction lower;

	public IntervalResultSparse(IntToDoubleFunction lower, IntToDoubleFunction upper)
	{
		this.upper = upper;
		this.lower = lower;
	}

	@Override public double getUpper(int state)
	{
		return upper.applyAsDouble(state);
	}

	@Override public double getLower(int state)
	{
		return lower.applyAsDouble(state);
	}

	@Override public double getAverage(int state)
	{
		return (upper.applyAsDouble(state) + lower.applyAsDouble(state)) / 2d;
	}
}
