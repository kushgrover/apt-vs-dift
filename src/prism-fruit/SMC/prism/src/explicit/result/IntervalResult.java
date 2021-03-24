package explicit.result;

import explicit.ModelCheckerResult;

import java.util.Arrays;

public interface IntervalResult
{
	double getUpper(int state);

	double getLower(int state);

	double getAverage(int state);

	static ModelCheckerResult toModelCheckerResult(ProcessingResult<IntervalResult> result, int numStates)
	{
		ModelCheckerResult res = new ModelCheckerResult();
		IntervalResult intervals = result.getResult();

		res.soln = new double[numStates];
		if (intervals instanceof IntervalResultSingle) {
			Arrays.fill(res.soln, ((IntervalResultSingle) intervals).getAverage());
		} else {
			Arrays.setAll(res.soln, intervals::getAverage);
		}
		res.numIters = result.getNumIters();
		res.timeTaken = result.getTimeTaken();
		return res;
	}
}
