package explicit.result;

public class GainBiasResultDense extends GainResultDense implements GainBiasResult
{
	private final double[] bias;

	public GainBiasResultDense(double[] gain, double[] bias)
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
