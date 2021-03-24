package simulator.method;

import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.RelOp;
import prism.PrismException;
import prism.PrismUtils;

/**
 * @author Przemysław Daca
 */
public class SPRTMethodError extends SPRTMethod
{
	// bound on the minimal transition probability
	private double minprob;
	// error: maximum probability that the sampler returns an arbitrary decision
	private double error;

	/**
	 * Use this constructor when the max. probability of finding false candidates is given.
	 */
	public SPRTMethodError(double alpha, double beta, double width, double falsecnd, double minprob)
	{
		super(alpha, beta, width);
		this.minprob = minprob;
		this.error = falsecnd;
	}

	@Override
	public SimulationMethod clone()
	{
		SPRTMethodError m = new SPRTMethodError(alpha, beta, delta, error, minprob);
		m.p0 = p0;
		m.p1 = p1;
		m.h0true = h0true;
		m.computedIterations = computedIterations;
		m.missingParameterComputed = missingParameterComputed;
		return m;
	}

	@Override
	public void setExpression(Expression expr) throws PrismException
	{
		super.setExpression(expr);
		RelOp relOp;
		if (expr instanceof ExpressionProb) {
			relOp = ((ExpressionProb) expr).getRelOp();
		}
		// For R properties...
		else if (expr instanceof ExpressionReward) {
			relOp = ((ExpressionReward) expr).getRelOp();
		}
		// Other (error)
		else {
			throw new PrismException("Cannot approximate " + expr + " using simulation");
		}

		/*
		if (relOp.isLowerBound()) {
			p0 = (1-error) * p0;
			p1 = (1-error) * p1 + error;
		}
		else {
			p0 = (1-error) * p0 + error;
			p1 = (1-error) * p1;
		}*/

		if (relOp.isLowerBound()) {
			p0 = p0 - error;
			p1 = p1 + error;
		} else {
			p0 = p0 + error;
			p1 = p1 - error;
		}

		if ((relOp.isLowerBound() && p0 <= p1) || (!relOp.isLowerBound() && p1 <= p0)) {
			throw new PrismException("Could not solve simulation parameters, indeference region to narrow!");
		}

	}

	@Override
	public String getParametersString()
	{
		if (!missingParameterComputed) {
			return "type I error=" + alpha + ", type II error=" + beta + ", width=" + delta + ", error=" + PrismUtils.formatDouble(3, error) + ", min. prob="
					+ minprob + ", iterations=unknown";
		} else {
			return "type I error=" + alpha + ", type II error=" + beta + ", width=" + delta + ", error=" + PrismUtils.formatDouble(3, error) + ", min. prob="
					+ minprob + ", iterations=" + computedIterations;
		}
	}

}
