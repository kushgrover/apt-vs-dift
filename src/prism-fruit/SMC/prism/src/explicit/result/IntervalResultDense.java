package explicit.result;

public class IntervalResultDense implements IntervalResult
{
	private final double[] upper;
	private final double[] lower;

	public IntervalResultDense(double[] upper, double[] lower)
	{
		this.upper = upper;
		this.lower = lower;
	}

	@Override public double getUpper(int state)
	{
		return upper[state];
	}

	@Override public double getLower(int state)
	{
		return lower[state];
	}

	@Override public double getAverage(int state)
	{
		return (upper[state] + lower[state]) / 2d;
	}

	public double[] getLower()
	{
		return lower;
	}

	public double[] getUpper()
	{
		return upper;
	}
}
