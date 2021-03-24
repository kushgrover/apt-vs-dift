package explicit.result;

import java.util.function.IntToDoubleFunction;

public class GainBiasResultSparse implements GainBiasResult
{
	private final IntToDoubleFunction gainProvider;
	private final IntToDoubleFunction biasProvider;

	public GainBiasResultSparse(IntToDoubleFunction gainProvider, IntToDoubleFunction biasProvider)
	{
		this.gainProvider = gainProvider;
		this.biasProvider = biasProvider;
	}

	@Override public double getGain(int state)
	{
		return gainProvider.applyAsDouble(state);
	}

	@Override public double getBias(int state)
	{
		return biasProvider.applyAsDouble(state);
	}
}
