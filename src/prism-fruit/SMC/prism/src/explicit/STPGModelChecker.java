//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Vojtech Forejt <vojtech.forejt@cs.ox.ac.uk> (University of Oxford)
//
//------------------------------------------------------------------------------
//
//	This file is part of PRISM.
//
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//==============================================================================

package explicit;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.rewards.STPGRewards;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import parser.ast.Expression;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import prism.PrismNotSupportedException;
import prism.PrismUtils;
import strat.MDStrategy;
import strat.MDStrategyArray;

import java.util.Arrays;
import java.util.Map;

/**
 * Explicit-state model checker for two-player stochastic games (STPGs).
 */
public class STPGModelChecker extends ProbModelChecker
{
	/**
	 * Create a new STPGModelChecker, inherit basic state from parent (unless null).
	 */
	public STPGModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	// Model checking functions

	@Override
	protected StateValues checkProbPathFormulaLTL(Model model, Expression expr, boolean qual, MinMax minMax, NatBitSet statesOfInterest) throws PrismException
	{
		throw new PrismNotSupportedException("LTL model checking not yet supported for stochastic games");
	}

	// Numerical computation functions

	/**
	 * Compute next=state probabilities.
	 * i.e. compute the probability of being in a state in {@code target} in the next step.
	 * @param stpg The STPG
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeNextProbs(STPG stpg, NatBitSet target, boolean min1, boolean min2) throws PrismException
	{
		ModelCheckerResult res = null;
		int n;
		double soln[], soln2[];
		long timer;

		timer = System.currentTimeMillis();

		// Store num states
		n = stpg.getNumStates();

		// Create/initialise solution vector(s)
		soln = Utils.intSetToDoubleArray(target, n);
		soln2 = new double[n];

		// Next-step probabilities
		stpg.mvMultMinMax(soln, min1, min2, soln2, null, null);

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln2;
		res.numIters = 1;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target}.
	 * @param stpg The STPG
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeReachProbs(STPG stpg, NatBitSet target, boolean min1, boolean min2) throws PrismException
	{
		return computeReachProbs(stpg, null, target, min1, min2, null, null);
	}

	/**
	 * Compute until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeUntilProbs(STPG stpg, NatBitSet remain, NatBitSet target, boolean min1, boolean min2) throws PrismException
	{
		return computeReachProbs(stpg, remain, target, min1, min2, null, null);
	}

	/**
	 * Compute reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (may be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	public ModelCheckerResult computeReachProbs(STPG stpg, NatBitSet remain, NatBitSet target, boolean min1, boolean min2, double init[], NatBitSet known)
			throws PrismException
	{
		ModelCheckerResult res = null;
		NatBitSet no, yes;
		int n, numYes, numNo;
		long timer, timerProb0, timerProb1;
		boolean genAdv;

		// Check for some unsupported combinations
		if (solnMethod == SolnMethod.VALUE_ITERATION && valIterDir == ValIterDir.ABOVE && !(precomp && prob0)) {
			throw new PrismException("Precomputation (Prob0) must be enabled for value iteration from above");
		}

		// Are we generating an optimal adversary?
		genAdv = exportAdv;

		// Start probabilistic reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting probabilistic reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		stpg.checkForDeadlocks(target);

		// Store num states
		n = stpg.getNumStates();

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null && !known.isEmpty()) {
			NatBitSet targetNew = target.clone();
			IntIterator iterator = known.iterator();
			while (iterator.hasNext()) {
				int i = iterator.nextInt();
				if (init[i] == 1.0) {
					targetNew.set(i);
				}
			}
			target = targetNew;
		}

		// Precomputation
		timerProb0 = System.currentTimeMillis();
		if (precomp && prob0) {
			no = prob0(stpg, remain, target, min1, min2);
		} else {
			no = NatBitSets.set();
		}
		timerProb0 = System.currentTimeMillis() - timerProb0;
		timerProb1 = System.currentTimeMillis();
		if (precomp && prob1 && !genAdv) {
			yes = prob1(stpg, remain, target, min1, min2);
		} else {
			yes = target.clone();
		}
		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numYes = yes.size();
		numNo = no.size();
		if (verbosity >= 1)
			mainLog.println("target=" + target.size() + ", yes=" + numYes + ", no=" + numNo + ", maybe=" + (n - (numYes + numNo)));

		// Compute probabilities
		switch (solnMethod) {
		case VALUE_ITERATION:
			res = computeReachProbsValIter(stpg, no, yes, min1, min2, init, known);
			break;
		case GAUSS_SEIDEL:
			res = computeReachProbsGaussSeidel(stpg, no, yes, min1, min2, init, known);
			break;
		default:
			throw new PrismException("Unknown STPG solution method " + solnMethod);
		}

		// Finished probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1)
			mainLog.println("Probabilistic reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		res.timeProb0 = timerProb0 / 1000.0;
		res.timePre = (timerProb0 + timerProb1) / 1000.0;

		return res;
	}

	/**
	 * Prob0 precomputation algorithm.
	 * i.e. determine the states of an STPG which, with min/max probability 0,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 * {@code min}=true gives Prob0E, {@code min}=false gives Prob0A.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public NatBitSet prob0(STPG stpg, NatBitSet remain, NatBitSet target, boolean min1, boolean min2)
	{
		int n, iters;
		NatBitSet u, soln, unknown;
		boolean u_done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting Prob0 (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Special case: no target states
		if (target.isEmpty()) {
			soln = NatBitSets.boundedFilledSet(stpg.getNumStates());
			return soln;
		}

		// Initialise vectors
		n = stpg.getNumStates();
		u = NatBitSets.boundedSet(n);
		soln = NatBitSets.boundedSet(n);

		// Determine set of states actually need to perform computation for
		unknown = NatBitSets.boundedFilledSet(n);
		unknown.andNot(target);
		if (remain != null)
			unknown.and(remain);

		// Fixed point loop
		iters = 0;
		u_done = false;
		// Least fixed point - should start from 0 but we optimise by
		// starting from 'target', thus bypassing first iteration
		u.or(target);
		soln.or(target);
		while (!u_done) {
			iters++;
			// Single step of Prob0
			stpg.prob0step(unknown, u, min1, min2, soln);
			// Check termination
			u_done = soln.equals(u);
			// u = soln
			u.clear();
			u.or(soln);
		}

		// Negate
		u.flip(0, n);

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Prob0 (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		return u;
	}

	/**
	 * Prob1 precomputation algorithm.
	 * i.e. determine the states of an STPG which, with min/max probability 1,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public NatBitSet prob1(STPG stpg, NatBitSet remain, NatBitSet target, boolean min1, boolean min2)
	{
		int n, iters;
		NatBitSet u, v, soln, unknown;
		boolean u_done, v_done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting Prob1 (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Special case: no target states
		if (target.isEmpty()) {
			return NatBitSets.boundedSet(stpg.getNumStates());
		}

		// Initialise vectors
		n = stpg.getNumStates();
		// Greatest fixed point
		u = NatBitSets.boundedFilledSet(n);
		v = NatBitSets.boundedSet(n);
		soln = NatBitSets.boundedSet(n);

		// Determine set of states actually need to perform computation for
		unknown = NatBitSets.boundedFilledSet(n);
		unknown.andNot(target);
		if (remain != null)
			unknown.and(remain);

		// Nested fixed point loop
		iters = 0;
		u_done = false;
		while (!u_done) {
			v_done = false;
			// Least fixed point - should start from 0 but we optimise by
			// starting from 'target', thus bypassing first iteration
			v.clear();
			v.or(target);
			soln.clear();
			soln.or(target);
			while (!v_done) {
				iters++;
				// Single step of Prob1
				stpg.prob1step(unknown, u, v, min1, min2, soln);
				// Check termination (inner)
				v_done = soln.equals(v);
				// v = soln
				v.clear();
				v.or(soln);
			}
			// Check termination (outer)
			u_done = v.equals(u);
			// u = v
			u.clear();
			u.or(v);
		}

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Prob1 (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		return u;
	}

	/**
	 * Compute reachability probabilities using value iteration.
	 * @param stpg The STPG
	 * @param no Probability 0 states
	 * @param yes Probability 1 states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (will be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachProbsValIter(STPG stpg, NatBitSet no, NatBitSet yes, boolean min1, boolean min2, double init[], NatBitSet known)
			throws PrismException
	{
		ModelCheckerResult res = null;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[], initVal;
		MDStrategy adv = null;
		boolean genAdv, done;
		long timer;

		// Are we generating an optimal adversary?
		genAdv = exportAdv;

		// Start value iteration
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector(s)
		soln = new double[n];
		soln2 = (init == null) ? new double[n] : init;

		// Initialise solution vectors. Use (where available) the following in order of preference:
		// (1) exact answer, if already known; (2) 1.0/0.0 if in yes/no; (3) passed in initial value; (4) initVal
		// where initVal is 0.0 or 1.0, depending on whether we converge from below/above.
		initVal = (valIterDir == ValIterDir.BELOW) ? 0.0 : 1.0;
		if (init != null) {
			if (known != null) {
				for (i = 0; i < n; i++)
					soln[i] = soln2[i] = known.contains(i) ? init[i] : yes.contains(i) ? 1.0 : no.contains(i) ? 0.0 : init[i];
			} else {
				for (i = 0; i < n; i++)
					soln[i] = soln2[i] = yes.contains(i) ? 1.0 : no.contains(i) ? 0.0 : init[i];
			}
		} else {
			for (i = 0; i < n; i++)
				soln[i] = soln2[i] = yes.contains(i) ? 1.0 : no.contains(i) ? 0.0 : initVal;
		}

		// Determine set of states actually need to compute values for
		unknown = NatBitSets.boundedFilledSet(n);
		unknown.andNot(yes);
		unknown.andNot(no);
		if (known != null)
			unknown.andNot(known);

		// Create/initialise adversary storage
		if (genAdv) {
			int[] advArray = new int[n];
			Arrays.fill(advArray, -1);
			adv = new MDStrategyArray(stpg, advArray);
		}

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			iters++;
			// Matrix-vector multiply and min/max ops
			stpg.mvMultMinMax(soln, min1, min2, soln2, unknown, genAdv ? adv : null);
			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Print adversary
		if (genAdv) {
			PrismLog out = new PrismFileLog(exportAdvFilename);
			for (i = 0; i < n; i++) {
				out.println(i + " " + (adv.isChoiceDefined(i) ? stpg.getAction(i, adv.getChoiceIndex(i)) : "-"));
			}
			out.println();
			out.close();
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute reachability probabilities using Gauss-Seidel.
	 * @param stpg The STPG
	 * @param no Probability 0 states
	 * @param yes Probability 1 states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (will be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachProbsGaussSeidel(STPG stpg, NatBitSet no, NatBitSet yes, boolean min1, boolean min2, double init[], NatBitSet known)
			throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], initVal, maxDiff;
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector
		soln = (init == null) ? new double[n] : init;

		// Initialise solution vector. Use (where available) the following in order of preference:
		// (1) exact answer, if already known; (2) 1.0/0.0 if in yes/no; (3) passed in initial value; (4) initVal
		// where initVal is 0.0 or 1.0, depending on whether we converge from below/above.
		initVal = (valIterDir == ValIterDir.BELOW) ? 0.0 : 1.0;
		if (init != null) {
			if (known != null) {
				for (i = 0; i < n; i++)
					soln[i] = known.contains(i) ? init[i] : yes.contains(i) ? 1.0 : no.contains(i) ? 0.0 : init[i];
			} else {
				for (i = 0; i < n; i++)
					soln[i] = yes.contains(i) ? 1.0 : no.contains(i) ? 0.0 : init[i];
			}
		} else {
			for (i = 0; i < n; i++)
				soln[i] = yes.contains(i) ? 1.0 : no.contains(i) ? 0.0 : initVal;
		}

		// Determine set of states actually need to compute values for
		unknown = NatBitSets.boundedFilledSet(n);
		unknown.andNot(yes);
		unknown.andNot(no);
		if (known != null)
			unknown.andNot(known);

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			iters++;
			// Matrix-vector multiply and min/max ops
			maxDiff = stpg.mvMultGSMinMax(soln, min1, min2, unknown, termCrit == TermCrit.ABSOLUTE);
			// Check termination
			done = maxDiff < termCritParam;
		}

		// Finished Gauss-Seidel
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Construct strategy information for min/max reachability probabilities.
	 * (More precisely, list of indices of player 1 choices resulting in min/max.)
	 * (Note: indices are guaranteed to be sorted in ascending order.)
	 * @param stpg The STPG
	 * @param state The state to generate strategy info for
	 * @param target The set of target states to reach
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param lastSoln Vector of probabilities from which to recompute in one iteration
	 */
	public IntList probReachStrategy(STPG stpg, int state, NatBitSet target, boolean min1, boolean min2, double lastSoln[]) throws PrismException
	{
		double val = stpg.mvMultMinMaxSingle(state, lastSoln, min1, min2);
		return stpg.mvMultMinMaxSingleChoices(state, lastSoln, min1, min2, val);
	}

	/**
	 * Compute bounded reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target} within k steps.
	 * @param stpg The STPG
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeBoundedReachProbs(STPG stpg, NatBitSet target, int k, boolean min1, boolean min2) throws PrismException
	{
		return computeBoundedReachProbs(stpg, null, target, k, min1, min2, null, null);
	}

	/**
	 * Compute bounded until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeBoundedUntilProbs(STPG stpg, NatBitSet remain, NatBitSet target, int k, boolean min1, boolean min2) throws PrismException
	{
		return computeBoundedReachProbs(stpg, remain, target, k, min1, min2, null, null);
	}

	/**
	 * Compute bounded reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Initial solution vector - pass null for default
	 * @param results Optional array of size k+1 to store (init state) results for each step (null if unused)
	 */
	public ModelCheckerResult computeBoundedReachProbs(STPG stpg, NatBitSet remain, NatBitSet target, int k, boolean min1, boolean min2, double init[],
			double results[]) throws PrismException
	{
		// TODO: implement until

		ModelCheckerResult res = null;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		long timer;

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting bounded probabilistic reachability...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector(s)
		soln = new double[n];
		soln2 = (init == null) ? new double[n] : init;

		// Initialise solution vectors. Use passed in initial vector, if present
		if (init != null) {
			for (i = 0; i < n; i++)
				soln[i] = soln2[i] = target.contains(i) ? 1.0 : init[i];
		} else {
			for (i = 0; i < n; i++)
				soln[i] = soln2[i] = target.contains(i) ? 1.0 : 0.0;
		}
		// Store intermediate results if required
		// (compute min/max value over initial states for first step)
		if (results != null) {
			results[0] = Utils.minMaxOverArraySubset(soln2, stpg.getInitialStates(), min2);
		}

		// Start iterations
		iters = 0;
		while (iters < k) {
			iters++;
			// Matrix-vector multiply and min/max ops
			stpg.mvMultMinMax(soln, min1, min2, soln2, NatBitSets.ensureBounded(target, stpg.getNumStates()).complement(), null);
			// Store intermediate results if required
			// (compute min/max value over initial states for this step)
			if (results != null) {
				results[iters] = Utils.minMaxOverArraySubset(soln2, stpg.getInitialStates(), min2);
			}
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Print vector (for debugging)
		//mainLog.println(soln);

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Bounded probabilistic reachability (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.lastSoln = soln2;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		res.timePre = 0.0;
		return res;
	}

	/**
	 * Compute expected reachability rewards.
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max rewards for player 1 (true=min, false=max)
	 * @param min2 Min or max rewards for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeReachRewards(STPG stpg, STPGRewards rewards, NatBitSet target, boolean min1, boolean min2) throws PrismException
	{
		return computeReachRewards(stpg, rewards, target, min1, min2, null, null);
	}

	/**
	 * Compute expected reachability rewards.
	 * i.e. compute the min/max reward accumulated to reach a state in {@code target}.
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max rewards for player 1 (true=min, false=max)
	 * @param min2 Min or max rewards for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (may be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	public ModelCheckerResult computeReachRewards(STPG stpg, STPGRewards rewards, NatBitSet target, boolean min1, boolean min2, double init[], NatBitSet known)
			throws PrismException
	{
		ModelCheckerResult res = null;
		NatBitSet inf;
		int n, numTarget, numInf;
		long timer, timerProb1;

		// Start expected reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting expected reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		stpg.checkForDeadlocks(target);

		// Store num states
		n = stpg.getNumStates();

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null && !known.isEmpty()) {
			NatBitSet targetNew = target.clone();
			IntIterator iterator = known.iterator();
			while (iterator.hasNext()) {
				int i = iterator.nextInt();
				if (init[i] == 1.0) {
					targetNew.set(i);
				}
			}
			target = targetNew;
		}

		// Precomputation (not optional)
		timerProb1 = System.currentTimeMillis();
		inf = prob1(stpg, null, target, !min1, !min2);
		inf.flip(0, n);
		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numTarget = target.size();
		numInf = inf.size();
		if (verbosity >= 1)
			mainLog.println("target=" + numTarget + ", inf=" + numInf + ", rest=" + (n - (numTarget + numInf)));

		// Compute rewards
		switch (solnMethod) {
		case VALUE_ITERATION:
			res = computeReachRewardsValIter(stpg, rewards, target, inf, min1, min2, init, known);
			break;
		default:
			throw new PrismException("Unknown STPG solution method " + solnMethod);
		}

		// Finished expected reachability
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1)
			mainLog.println("Expected reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		res.timePre = timerProb1 / 1000.0;

		return res;
	}

	/**
	 * Compute expected reachability rewards using value iteration.
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param inf States for which reward is infinite
	 * @param min1 Min or max rewards for player 1 (true=min, false=max)
	 * @param min2 Min or max rewards for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (will be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachRewardsValIter(STPG stpg, STPGRewards rewards, NatBitSet target, NatBitSet inf, boolean min1, boolean min2,
			double init[], NatBitSet known) throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector(s)
		soln = new double[n];
		soln2 = (init == null) ? new double[n] : init;

		// Initialise solution vectors. Use (where available) the following in order of preference:
		// (1) exact answer, if already known; (2) 0.0/infinity if in target/inf; (3) passed in initial value; (4) 0.0
		if (init != null) {
			if (known != null) {
				for (i = 0; i < n; i++)
					soln[i] = soln2[i] = known.contains(i) ? init[i] : target.contains(i) ? 0.0 : inf.contains(i) ? Double.POSITIVE_INFINITY : init[i];
			} else {
				for (i = 0; i < n; i++)
					soln[i] = soln2[i] = target.contains(i) ? 0.0 : inf.contains(i) ? Double.POSITIVE_INFINITY : init[i];
			}
		} else {
			for (i = 0; i < n; i++)
				soln[i] = soln2[i] = target.contains(i) ? 0.0 : inf.contains(i) ? Double.POSITIVE_INFINITY : 0.0;
		}

		// Determine set of states actually need to compute values for
		unknown = NatBitSets.boundedFilledSet(n);
		unknown.andNot(target);
		unknown.andNot(inf);
		if (known != null)
			unknown.andNot(known);

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			//mainLog.println(soln);
			iters++;
			// Matrix-vector multiply and min/max ops
			stpg.mvMultRewMinMax(soln, rewards, min1, min2, soln2, unknown, null);
			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Simple test program.
	 */
	public static void main(String args[])
	{
		STPGModelChecker mc;
		STPGAbstrSimple stpg;
		ModelCheckerResult res;
		NatBitSet target;
		Map<String, NatBitSet> labels;
		boolean min1 = true, min2 = true;
		try {
			mc = new STPGModelChecker(null);
			stpg = new STPGAbstrSimple();
			stpg.buildFromPrismExplicit(args[0]);
			stpg.addInitialState(0);
			//System.out.println(stpg);
			labels = StateModelChecker.loadLabelsFile(args[1]);
			//System.out.println(labels);
			target = labels.get(args[2]);
			if (target == null)
				throw new PrismException("Unknown label \"" + args[2] + "\"");
			for (int i = 3; i < args.length; i++) {
				if (args[i].equals("-minmin")) {
					min1 = true;
					min2 = true;
				} else if (args[i].equals("-maxmin")) {
					min1 = false;
					min2 = true;
				} else if (args[i].equals("-minmax")) {
					min1 = true;
					min2 = false;
				} else if (args[i].equals("-maxmax")) {
					min1 = false;
					min2 = false;
				}
			}
			//stpg.exportToDotFile("stpg.dot", target);
			//stpg.exportToPrismExplicit("stpg");
			res = mc.computeReachProbs(stpg, target, min1, min2);
			System.out.println(res.soln[0]);
		} catch (PrismException e) {
			System.out.println(e);
		}
	}
}
