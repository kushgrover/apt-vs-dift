package simulator.method;

import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.RelOp;
import prism.PrismException;
import simulator.sampler.Sampler;
import simulator.sampler.SamplerCandidateLTL;
import simulator.sampler.SamplerDouble;

/**
 * Written by Maxi on 14.10.20 as a quick hack for the TACAS paper with Stefan Kiefer, Javier and Jan.
 */

public class MonitorMethod extends SimulationMethod
{
	//Variant is used to determine whether we are cautious or bold monitor
	double variant;

	public MonitorMethod(double variant)
	{
		this.variant = variant;
	}

	@Override
	public String getName()
	{
		return (variant > 0.5) ? "New Bold (sophisticated) Monitor" : "New Cautious (naive) Monitor";
	}

	@Override
	public String getFullName()
	{
		return getName();
	}

	@Override
	public void reset()
	{
		//nothing to do here
	}

	@Override
	public void computeMissingParameterBeforeSim() throws PrismException
	{
		// Nothing to do
	}

	/**
	 * Copy pasted from SPRTMethod, hope that works out; but we are doing sth that is a subset of SPRT, so should be fine.
	 */
	@Override
	public void setExpression(Expression expr) throws PrismException
	{
		Expression bound;
		RelOp relOp;
		double theta;

		// For P properties...
		if (expr instanceof ExpressionProb) {
			bound = ((ExpressionProb) expr).getProb();
			relOp = ((ExpressionProb) expr).getRelOp();
		}
		// For R properties...
		else if (expr instanceof ExpressionReward) {
			bound = ((ExpressionReward) expr).getReward();
			relOp = ((ExpressionReward) expr).getRelOp();
		}
		// Other (error)
		else {
			throw new PrismException("Cannot approximate " + expr + " using simulation");
		}

//		// There must be a probability/reward bound to use SPRT
		// Maxi: But for monitor, we don't return a model checking result anyway, so we probably don't care.
//		if (bound == null) {
//			throw new PrismException("Cannot use SPRT on a quantitative (=?) property");
//		} else {
//			theta = bound.evaluateDouble();
//			// Check of out-of-range errors
//			if (theta - delta <= 0) {
//				String s = "Indifference for SPRT method (" + delta + ") is too wide";
//				s += " (bound " + theta + " is too close to 0)";
//				throw new PrismException(s);
//			}
//			if (theta + delta > 1 && expr instanceof ExpressionProb) {
//				String s = "Indifference for SPRT method (" + delta + ") is too wide";
//				s += " (bound " + theta + " is too close to 1)";
//				throw new PrismException(s);
//			}
//			// Set p0/p1 values for two hypotheses
//			// Default case is that H0 means probability/reward >= theta + delta
//			if (relOp.isLowerBound()) {
//				p0 = theta + delta;
//				p1 = theta - delta;
//			}
//			// If the bound is reversed, just swap p0/p1
//			else {
//				p0 = theta - delta;
//				p1 = theta + delta;
//			}
//		}
	}

	@Override
	public void computeMissingParameterAfterSim()
	{
		// Nothing to do (this is done in shouldStopNow)
	}

	@Override
	public Object getMissingParameter() throws PrismException
	{
		return "getMissingParameter in MonitorMethod doesn't do anything.";
	}

	@Override
	public String getParametersString()
	{
		return "Variant: " + variant + ", hence we are a " + getName();
	}

	@Override
	public boolean shouldStopNow(int iters, Sampler sampler)
	{
		if (!(sampler instanceof SamplerCandidateLTL)){
			System.err.println("MonitorMethod needs to be run with SamplerCandidateLTL");
			System.exit(-1);
		}
		if (sampler.isCurrentValueKnown()){
			return (boolean) sampler.getCurrentValue();
		}
		else{
			return false;
		}
	}

	@Override
	public int getProgress(int iters, Sampler sampler)
	{
		// No good measure of progress unfortunately
		return 0;
	}

	@Override
	public Object getResult(Sampler sampler) throws PrismException
	{
		//TODO: Maybe change?
		return true;
	}

	@Override
	public String getResultExplanation(Sampler sampler){
		//This monitor only works for LTL, so we can "safely" cast.
		SamplerCandidateLTL sam = ((SamplerCandidateLTL) sampler);
		//get the number of resets. Since last simulation succeeded (else we wouldn't have stopped), this is NumSamples-1
		int resets = sam.getNumSamples() - 1;
		//get the total steps the simulator took minus the steps of the last run, hence the number of steps in all reset runs
		long stepsInAllResetRuns = sam.getCanTotal() + sam.getTransTotal() - sam.getLastRunSteps();
		return "Resets: " + resets + " ;stepsInResetRuns: " + stepsInAllResetRuns;
	}

	@Override
	public SimulationMethod clone()
	{
		return new MonitorMethod(variant);
	}



}
