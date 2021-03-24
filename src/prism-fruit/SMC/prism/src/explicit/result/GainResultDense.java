package explicit.result;

public class GainResultDense implements GainResult
{
	private final double[] gain;

	public GainResultDense(double[] gain)
	{
		this.gain = gain;
	}

	@Override public double getGain(int state)
	{
		return gain[state];
	}

	public double[] getGain()
	{
		return gain;
	}
}
