package explicit.result;

import explicit.ModelCheckerResult;

import java.util.Arrays;

public interface GainResult
{
	static ModelCheckerResult toModelCheckerResult(ProcessingResult<? extends GainResult> result, int numStates)
	{
		ModelCheckerResult res = new ModelCheckerResult();
		GainResult gain = result.getResult();

		if (gain instanceof GainResultDense) {
			res.soln = ((GainResultDense) gain).getGain();
		} else {
			res.soln = new double[numStates];
			if (gain instanceof GainResultSingle) {
				Arrays.fill(res.soln, ((GainResultSingle) gain).getGain());
			} else {
				Arrays.setAll(res.soln, gain::getGain);
			}
		}
		res.numIters = result.getNumIters();
		res.timeTaken = result.getTimeTaken();
		return res;
	}

	double getGain(int state);
}
