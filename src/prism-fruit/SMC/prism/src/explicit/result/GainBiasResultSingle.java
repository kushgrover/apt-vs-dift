package explicit.result;

public class GainBiasResultSingle extends GainResultSingle implements GainBiasResult
{
	private final double[] bias;

	public GainBiasResultSingle(double gain, double[] bias)
	{
		super(gain);
		this.bias = bias;
	}

	@Override public double getBias(int state)
	{
		return bias[state];
	}

	public double[] getBias()
	{
		return bias;
	}
}
