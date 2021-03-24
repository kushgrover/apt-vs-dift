package simulator.method;

import cern.colt.Timer;
import org.apache.commons.math3.distribution.BinomialDistribution;
import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.RelOp;
import prism.Pair;
import prism.PrismException;
import simulator.sampler.Sampler;

/**
 * Simulation method for single sampling plan.
 *
 * @author Przemysław Daca
 */
public class SSPMethod extends SimulationMethod
{
	// Hypothesis H0: probability/expectation >= theta + delta
	// Hypothesis H1: probability/expectation <= theta - delta

	// TODO make this some parameter
	private static int maxSampleSize = 1000000;
	// SPRT parameters
	// alpha: probability of type I error (false acceptance of H1)
	protected double alpha;
	// beta: probability of type II error (false acceptance of H0)
	protected double beta;
	// width: half the width of the indifference region
	protected double width;
	// Probability/expectation for hypothesis 0/1 (used for likelihood ratio)
	protected double p0;
	protected double p1;
	// number of iterations and the critical value
	protected int sampleSize;
	protected int criticalVal;
	// does H0 represent that the property holds?
	protected boolean propertyIsH0;
	Timer paramTmr;
	// decision made
	private boolean h0true;

	public SSPMethod(double alpha, double beta, double width)
	{
		this.alpha = alpha;
		this.beta = beta;
		this.width = width;
		this.p0 = 0;
		this.p1 = 0;
		this.sampleSize = 0;
		this.criticalVal = 0;
		this.paramTmr = new Timer();
	}

	@Override
	public String getName()
	{
		return "SSP";
	}

	@Override
	public String getFullName()
	{
		return "Single Sampling Plan";
	}

	@Override
	public void reset()
	{
		h0true = false;
	}

	@Override
	public void computeMissingParameterBeforeSim() throws PrismException
	{
	}

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

		// There must be a probability/reward bound
		if (bound == null) {
			throw new PrismException("Cannot use SPRT on a quantitative (=?) property");
		} else {
			theta = bound.evaluateDouble();
			// Check of out-of-range errors
			if (theta - width <= 0) {
				String s = "Indifference for SSP method (" + width + ") is too wide";
				s += " (bound " + theta + " is too close to 0)";
				throw new PrismException(s);
			}
			if (theta + width > 1 && expr instanceof ExpressionProb) {
				String s = "Indifference for SSP method (" + width + ") is too wide";
				s += " (bound " + theta + " is too close to 1)";
				throw new PrismException(s);
			}

			p0 = theta + width;
			p1 = theta - width;

			if (relOp.isLowerBound()) {
				propertyIsH0 = true;
			} else {
				// property is P[phi]<= p_1 instead of P[phi]>= p_0
				propertyIsH0 = false;
				double tmp = alpha;
				alpha = beta;
				beta = tmp;
			}
		}

		computeParameters();
	}

	/**
	 * Calculates the samples size and the critical value.
	 *
	 * @throws PrismException
	 */
	protected void computeParameters() throws PrismException
	{

		paramTmr.start();
		Pair<Integer, Integer> res = binarySearch(alpha, beta, p0, p1);
		if (res == null) {
			paramTmr.stop();
			throw new PrismException("Couldn't compute simulation parameters, for samples size <" + maxSampleSize);
		} else {
			sampleSize = res.first;
			criticalVal = res.second;
		}
		paramTmr.stop();
	}

	/**
	 * Find the sample size and critical value by binary search algorithm.
	 */
	private Pair<Integer, Integer> binarySearch(double alpha, double beta, double p0, double p1)
	{

		// do binary search on n;
		int nmin = 1;
		int nmax = maxSampleSize;
		int n, c;

		do {
			n = (nmax - nmin) / 2 + nmin;

			c = solveCV(n, p0, p1);

			if (c == -1) {
				nmin = n;
			} else {
				nmax = n;
			}

		} while ((nmax - nmin > 1));

		c = solveCV(nmax, p0, p1);

		if (c != -1) {
			return new Pair<>(nmax, c);
		} else {
			return null;
		}
	}

	/**
	 * Finds the critical value for the given parameters. Returns -1 if no solution exists.
	 */
	private int solveCV(int n, double p0, double p1)
	{
		// find the lowest value of c that the satisfies Type-I bound
		BinomialDistribution dlow = new BinomialDistribution(n, p0);
		int c = dlow.inverseCumulativeProbability(alpha) - 1;

		// check if c satisfies Type-II bound
		BinomialDistribution dhigh = new BinomialDistribution(n, p1);
		double prob = dhigh.cumulativeProbability(c);
		if (prob > 1 - beta) {
			return c;
		} else {
			return -1;
		}
	}

	@Override
	public void computeMissingParameterAfterSim()
	{
	}

	@Override
	public Object getMissingParameter() throws PrismException
	{
		return sampleSize;
	}

	@Override
	public String getParametersString()
	{
		if (sampleSize <= 0) {
			return "type I error=" + alpha + ", type II error=" + beta + ", width=" + width + ", samples=?, critical value=?";
		} else {
			return "type I error=" + alpha + ", type II error=" + beta + ", width=" + width + ", samples=" + sampleSize + ", critical value=" + criticalVal;
		}
	}

	@Override
	public boolean shouldStopNow(int iters, Sampler sampler)
	{
		if (iters == sampleSize) {
			// compute the decision
			double sat = sampler.getMeanValue() * sampleSize;
			h0true = (sat >= criticalVal);
		}

		return (iters > sampleSize);
	}

	@Override
	public int getProgress(int iters, Sampler sampler)
	{
		return 10 * ((int) (10 * iters) / sampleSize);

	}

	@Override
	public Object getResult(Sampler sampler) throws PrismException
	{
		return h0true ^ !propertyIsH0;
	}

	@Override
	public String getResultExplanation(Sampler sampler) throws PrismException
	{
		return "time to compute simulation params " + paramTmr.toString();
	}

	@Override

	public SimulationMethod clone()
	{
		SSPMethod clone = new SSPMethod(alpha, beta, width);
		clone.p0 = p0;
		clone.p1 = p1;
		clone.sampleSize = sampleSize;
		clone.criticalVal = criticalVal;
		clone.propertyIsH0 = propertyIsH0;

		return clone;
	}

}
