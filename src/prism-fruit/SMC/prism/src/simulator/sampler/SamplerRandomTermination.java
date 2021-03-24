package simulator.sampler;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.doubles.DoubleListIterator;
import parser.State;
import parser.ast.Expression;
import parser.ast.ExpressionTemporal;
import prism.PrismException;
import prism.PrismLangException;
import simulator.Path;
import simulator.TransitionList;

import java.util.Random;

/**
 * Sampler for unbounded until which terminates randomly. Shorter paths that satisfy the property have more weight.
 * Based on "Statistical Verification of Probabilistic Properties with Unbounded Until" Younes'11.
 * <p>
 * Note: although this sampler doesn't assign 0-1 values to paths, but implementation-wise it's easier to think of this sampler as boolean.
 *
 * @author Przemysław Daca
 */
public class SamplerRandomTermination extends SamplerBoolean
{
	private final Random rnd;
	// Value of current path
	protected boolean value;
	// Whether the actual value should be the negation of 'value'
	protected boolean negated = false;
	// Stats over all paths
	protected int numSamples;
	protected double sumValues;        // x_1 + ... + x_n
	private Expression left;
	private Expression right;
	private double pterm;
	private long pathsize;
	private DoubleList values;   // x_1, ..., x_n

	public SamplerRandomTermination(ExpressionTemporal expr, double probTerm) throws PrismException
	{
		this.pterm = probTerm;
		this.rnd = new Random(System.currentTimeMillis());
		// Make sure expression is of the correct type
		// Then extract other required info
		if (expr.getOperator() != ExpressionTemporal.P_U) {
			throw new PrismException("Error creating Sampler");
		}
		left = expr.getOperand1();
		right = expr.getOperand2();

		values = new DoubleArrayList();
		// Initialise sampler info
		reset();
		resetStats();
	}

	@Override
	public void reset()
	{
		valueKnown = false;
		value = false;

	}

	@Override
	public void resetStats()
	{
		numSamples = 0;
		sumValues = 0;
		pathsize = 0;
		values = new DoubleArrayList();
	}

	@Override
	public boolean update(Path path, TransitionList transList) throws PrismLangException
	{

		pathsize = path.size();

		// If the answer is already known we should do nothing
		if (valueKnown) {
			return true;
		}

		State currentState = path.getCurrentState();
		// Have we reached the target (i.e. RHS of until)?
		if (right.evaluateBoolean(currentState)) {
			valueKnown = true;
			value = true;
		}
		// Or, if not, have we violated the LHS of the until?
		else if (!left.evaluateBoolean(currentState)) {
			valueKnown = true;
			value = false;
		}

		// Przemek: I delete deadlock detection, since it's not really statistical model checking
		// Or, if we are now at a deadlock/self-loop
		/*else if (transList != null && (transList.isDeadlock() || transList.isDeterministicSelfLoop(currentState))) {
			valueKnown = true;
			value = false;
		}*/

		if (!valueKnown) {
			// terminate with given probability;
			if (rnd.nextDouble() < pterm) {
				valueKnown = true;
				value = false;
			}
		}

		return valueKnown;
	}

	@Override
	public void updateStats()
	{
		numSamples++;
		// XOR: value && !negated || !value && negated
		double xi = 0;

		if (value != negated) {
			xi = Math.pow(1 - pterm, -pathsize);
			sumValues += xi;
		}

		values.add(xi);
	}

	@Override
	public Object getCurrentValue()
	{
		return value != negated;
	}

	@Override
	public double getMeanValue()
	{
		return sumValues / numSamples;
	}

	@Override
	public double getVariance()
	{
		if (numSamples < 1) {
			return 0.0;
		}

		double sumOfSqr = 0;
		double mean = getMeanValue();

		DoubleListIterator iterator = values.iterator();
		while (iterator.hasNext()) {
			double xi = iterator.nextIndex();
			double sqr = (xi - mean);
			sumOfSqr += sqr * sqr;
		}

		double variance = (1 + sumOfSqr) / numSamples;
		return variance;
	}

	@Override
	public double getLikelihoodRatio(double p1, double p0) throws PrismException
	{
		throw new PrismException("LR not implemented");
	}
}
