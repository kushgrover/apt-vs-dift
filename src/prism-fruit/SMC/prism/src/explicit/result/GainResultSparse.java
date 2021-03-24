package explicit.result;

import java.util.function.IntToDoubleFunction;

public class GainResultSparse implements GainResult
{
	private final IntToDoubleFunction gainProvider;

	public GainResultSparse(IntToDoubleFunction gainProvider)
	{
		this.gainProvider = gainProvider;
	}

	@Override public double getGain(int state)
	{
		return gainProvider.applyAsDouble(state);
	}
}
