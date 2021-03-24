package explicit.result;

public class IntervalResultSingle implements IntervalResult
{
	private final double upper;
	private final double lower;

	public IntervalResultSingle(double lower, double upper)
	{
		this.upper = upper;
		this.lower = lower;
	}

	@Override public double getUpper(int state)
	{
		return upper;
	}

	@Override public double getLower(int state)
	{
		return lower;
	}

	@Override public double getAverage(int state)
	{
		return (upper + lower) / 2d;
	}

	public double getLower()
	{
		return lower;
	}

	public double getUpper()
	{
		return upper;
	}

	public double getAverage()
	{
		return (upper + lower) / 2d;
	}
}
