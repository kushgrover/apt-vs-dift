package simulator.method;

import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.RelOp;
import prism.PrismException;
import prism.PrismUtils;

/**
 * Single sampling plan with error correction.
 *
 * @author Przemysław Daca
 */
public class SSPMethodError extends SSPMethod
{
	// Hypothesis H0: probability/expectation >= theta + delta
	// Hypothesis H1: probability/expectation <= theta - delta

	// bound on the minimal transition probability
	private double minprob;
	// error: maximum probability that the sampler returns an arbitrary decision
	private double error;

	/**
	 * Use this constructor when the max. probability of finding false candidates is given.
	 */
	public SSPMethodError(double alpha, double beta, double width, double falsecnd, double minprob)
	{
		super(alpha, beta, width);
		this.minprob = minprob;
		this.error = falsecnd;
	}

	@Override
	public SimulationMethod clone()
	{
		SSPMethodError clone = new SSPMethodError(alpha, beta, width, error, minprob);
		clone.p0 = p0;
		clone.p1 = p1;
		clone.sampleSize = sampleSize;
		clone.criticalVal = criticalVal;
		clone.propertyIsH0 = propertyIsH0;

		return clone;
	}

	@Override
	public void setExpression(Expression expr) throws PrismException
	{
		Expression bound;
		RelOp relOp;
		double theta;

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

		/*p0 = (1-error) * p0;
		p1 = (1-error) * p1 + error;*/
		p0 = p0 - error;
		p1 = p1 + error;

		if (p0 <= p1) {
			throw new PrismException("Error: Could not solve simulation parameters!");
		}

		computeParameters();
	}

	@Override
	public String getParametersString()
	{
		if (sampleSize <= 0) {
			return "type I error=" + alpha + ", type II error=" + beta + ", width=" + width + ", error=" + PrismUtils.formatDouble(3, error) + ", min. prob="
					+ minprob + ", iterations=?, critical value=?";
		} else {
			return "type I error=" + alpha + ", type II error=" + beta + ", width=" + width + ", error=" + PrismUtils.formatDouble(3, error) + ", min. prob="
					+ minprob + ", iterations=" + sampleSize + ", critical value=" + criticalVal;
		}
	}
}
