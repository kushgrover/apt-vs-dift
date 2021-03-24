package explicit.result;

public class GainResultSingle implements GainResult
{
	private final double gain;

	public GainResultSingle(double gain)
	{
		this.gain = gain;
	}

	@Override public double getGain(int state)
	{
		return gain;
	}

	public double getGain()
	{
		return gain;
	}
}
