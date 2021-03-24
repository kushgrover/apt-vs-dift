//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import acceptance.AcceptanceReach;
import acceptance.AcceptanceType;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;
import common.FastUtils;
import common.Time;
import de.tum.in.naturals.Indices;
import de.tum.in.naturals.map.Nat2DoubleDenseArrayMap;
import de.tum.in.naturals.set.BoundedNatBitSet;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.result.GainBiasResult;
import explicit.result.GainBiasResultDense;
import explicit.result.GainBiasResultSingle;
import explicit.result.GainBiasResultSparse;
import explicit.result.GainResult;
import explicit.result.IntervalResult;
import explicit.result.IntervalResultDense;
import explicit.result.IntervalResultSingle;
import explicit.result.IntervalResultSparse;
import explicit.result.ProcessingResult;
import explicit.rewards.MCRewards;
import explicit.rewards.Rewards;
import it.unimi.dsi.fastutil.doubles.DoubleIterator;
import it.unimi.dsi.fastutil.doubles.DoubleLinkedOpenHashSet;
import it.unimi.dsi.fastutil.doubles.DoubleSet;
import it.unimi.dsi.fastutil.ints.Int2DoubleFunction;
import it.unimi.dsi.fastutil.ints.Int2DoubleLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntListIterator;
import org.ojalgo.function.aggregator.Aggregator;
import org.ojalgo.matrix.decomposition.LU;
import org.ojalgo.matrix.decomposition.MatrixDecomposition;
import org.ojalgo.matrix.decomposition.SingularValue;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.PrimitiveDenseStore;
import org.ojalgo.matrix.store.SparseStore;
import org.ojalgo.matrix.store.SparseStore2;
import org.ojalgo.matrix.task.SolverTask;
import org.ojalgo.matrix.task.TaskException;
import parser.VarList;
import parser.ast.Declaration;
import parser.ast.DeclarationIntUnbounded;
import parser.ast.Expression;
import prism.Prism;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import prism.PrismNotSupportedException;
import prism.PrismSettings;
import prism.PrismUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.function.DoubleToIntFunction;
import java.util.function.IntConsumer;
import java.util.function.IntPredicate;
import java.util.function.IntToDoubleFunction;
import java.util.function.IntUnaryOperator;
import java.util.function.Supplier;

import static com.google.common.base.Preconditions.checkArgument;
import static explicit.DTMC.Type.STRONGLY_UNICHAIN;
import static explicit.DTMC.Type.UNICHAIN;

/**
 * Explicit-state model checker for discrete-time Markov chains (DTMCs).
 */
public class DTMCModelChecker extends ProbModelChecker
{
	/**
	 * Create a new DTMCModelChecker, inherit basic state from parent (unless null).
	 */
	public DTMCModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	/**
	 * Simple test program.
	 */
	public static void main(String args[])
	{
		DTMCModelChecker mc;
		DTMCSimple dtmc;
		ModelCheckerResult res;
		try {
			// Two examples of building and solving a DTMC

			int version = 2;
			if (version == 1) {

				// 1. Read in from .tra and .lab files
				//    Run as: PRISM_MAINCLASS=explicit.DTMCModelChecker bin/prism dtmc.tra dtmc.lab target_label
				mc = new DTMCModelChecker(null);
				dtmc = new DTMCSimple();
				dtmc.buildFromPrismExplicit(args[0]);
				dtmc.addInitialState(0);
				//System.out.println(dtmc);
				Map<String, NatBitSet> labels = StateModelChecker.loadLabelsFile(args[1]);
				//System.out.println(labels);
				NatBitSet target = labels.get(args[2]);
				if (target == null)
					throw new PrismException("Unknown label \"" + args[2] + "\"");
				for (int i = 3; i < args.length; i++) {
					if (args[i].equals("-nopre"))
						mc.setPrecomp(false);
				}
				res = mc.computeReachProbs(dtmc, target);
				System.out.println(res.soln[0]);

			} else {

				// 2. Build DTMC directly
				//    Run as: PRISM_MAINCLASS=explicit.DTMCModelChecker bin/prism
				//    (example taken from p.14 of Lec 5 of http://www.prismmodelchecker.org/lectures/pmc/)
				mc = new DTMCModelChecker(null);
				dtmc = new DTMCSimple(6);
				dtmc.setProbability(0, 1, 0.1);
				dtmc.setProbability(0, 2, 0.9);
				dtmc.setProbability(1, 0, 0.4);
				dtmc.setProbability(1, 3, 0.6);
				dtmc.setProbability(2, 2, 0.1);
				dtmc.setProbability(2, 3, 0.1);
				dtmc.setProbability(2, 4, 0.5);
				dtmc.setProbability(2, 5, 0.3);
				dtmc.setProbability(3, 3, 1.0);
				dtmc.setProbability(4, 4, 1.0);
				dtmc.setProbability(5, 5, 0.3);
				dtmc.setProbability(5, 4, 0.7);
				System.out.println(dtmc);
				NatBitSet target = NatBitSets.set();
				target.set(4);
				NatBitSet remain = NatBitSets.set();
				remain.set(1);
				remain.flip(0, 6);
				System.out.println(target);
				System.out.println(remain);
				res = mc.computeUntilProbs(dtmc, remain, target);
				System.out.println(res.soln[0]);
			}
		} catch (PrismException e) {
			System.out.println(e);
		}
	}

	// Model checking functions

	@Override
	protected StateValues checkProbPathFormulaLTL(Model model, Expression expr, boolean qual, MinMax minMax, NatBitSet statesOfInterest) throws PrismException
	{
		LTLModelChecker mcLtl;
		StateValues probsProduct, probs;
		LTLModelChecker.LTLProduct<DTMC> product;
		DTMCModelChecker mcProduct;

		// For LTL model checking routines
		mcLtl = new LTLModelChecker(this);

		// Build product of Markov chain and automaton
		AcceptanceType[] allowedAcceptance = {
				AcceptanceType.RABIN,
				AcceptanceType.REACH,
				AcceptanceType.GENERIC
		};
		product = mcLtl.constructProductMC(this, (DTMC) model, expr, statesOfInterest, allowedAcceptance);

		// Output product, if required
		if (getExportProductTrans()) {
			mainLog.println("\nExporting product transition matrix to file \"" + getExportProductTransFilename() + "\"...");
			product.getProductModel().exportToPrismExplicitTra(getExportProductTransFilename());
		}
		if (getExportProductStates()) {
			mainLog.println("\nExporting product state space to file \"" + getExportProductStatesFilename() + "\"...");
			PrismFileLog out = new PrismFileLog(getExportProductStatesFilename());
			VarList newVarList = (VarList) modulesFile.createVarList().clone();
			String daVar = "_da";
			while (newVarList.getIndex(daVar) != -1) {
				daVar = "_" + daVar;
			}
			newVarList.addVar(0, new Declaration(daVar, new DeclarationIntUnbounded()), 1, null);
			product.getProductModel().exportStates(Prism.EXPORT_PLAIN, newVarList, out);
			out.close();
		}

		// Find accepting states + compute reachability probabilities
		NatBitSet acc;
		if (product.getAcceptance() instanceof AcceptanceReach) {
			mainLog.println("\nSkipping BSCC computation since acceptance is defined via goal states...");
			acc = ((AcceptanceReach) product.getAcceptance()).getGoalStates();
		} else {
			mainLog.println("\nFinding accepting BSCCs...");
			acc = mcLtl.findAcceptingBSCCs(product.getProductModel(), product.getAcceptance());
		}
		mainLog.println("\nComputing reachability probabilities...");
		mcProduct = new DTMCModelChecker(this);
		mcProduct.inheritSettings(this);
		ModelCheckerResult res = mcProduct.computeReachProbs(product.getProductModel(), acc);
		probsProduct = StateValues.createFromDoubleArray(res.soln, product.getProductModel());

		// Output vector over product, if required
		if (getExportProductVector()) {
			mainLog.println("\nExporting product solution vector matrix to file \"" + getExportProductVectorFilename() + "\"...");
			PrismFileLog out = new PrismFileLog(getExportProductVectorFilename());
			probsProduct.print(out, false, false, false, false);
			out.close();
		}

		// Mapping probabilities in the original model
		probs = product.projectToOriginalModel(probsProduct);
		probsProduct.clear();

		return probs;
	}

	/**
	 * Compute rewards for a co-safe LTL reward operator.
	 */
	@Override protected StateValues checkRewardCoSafeLTL(Model model, Rewards modelRewards, Expression expr, MinMax minMax, NatBitSet statesOfInterest)
			throws PrismException
	{
		LTLModelChecker mcLtl;
		MCRewards productRewards;
		StateValues rewardsProduct, rewards;
		DTMCModelChecker mcProduct;
		LTLModelChecker.LTLProduct<DTMC> product;

		// For LTL model checking routines
		mcLtl = new LTLModelChecker(this);

		// Build product of Markov chain and automaton
		AcceptanceType[] allowedAcceptance = {
				AcceptanceType.RABIN,
				AcceptanceType.REACH
		};
		product = mcLtl.constructProductMC(this, (DTMC) model, expr, statesOfInterest, allowedAcceptance);

		// Adapt reward info to product model
		productRewards = ((MCRewards) modelRewards).liftFromModel(product);

		// Output product, if required
		if (getExportProductTrans()) {
			mainLog.println("\nExporting product transition matrix to file \"" + getExportProductTransFilename() + "\"...");
			product.getProductModel().exportToPrismExplicitTra(getExportProductTransFilename());
		}
		if (getExportProductStates()) {
			mainLog.println("\nExporting product state space to file \"" + getExportProductStatesFilename() + "\"...");
			PrismFileLog out = new PrismFileLog(getExportProductStatesFilename());
			VarList newVarList = (VarList) modulesFile.createVarList().clone();
			String daVar = "_da";
			while (newVarList.getIndex(daVar) != -1) {
				daVar = "_" + daVar;
			}
			newVarList.addVar(0, new Declaration(daVar, new DeclarationIntUnbounded()), 1, null);
			product.getProductModel().exportStates(Prism.EXPORT_PLAIN, newVarList, out);
			out.close();
		}

		// Find accepting states + compute reachability rewards
		NatBitSet acc;
		if (product.getAcceptance() instanceof AcceptanceReach) {
			mainLog.println("\nSkipping BSCC computation since acceptance is defined via goal states...");
			acc = ((AcceptanceReach) product.getAcceptance()).getGoalStates();
		} else {
			mainLog.println("\nFinding accepting BSCCs...");
			acc = mcLtl.findAcceptingBSCCs(product.getProductModel(), product.getAcceptance());
		}
		mainLog.println("\nComputing reachability probabilities...");
		mcProduct = new DTMCModelChecker(this);
		mcProduct.inheritSettings(this);
		ModelCheckerResult res = mcProduct.computeReachRewards(product.getProductModel(), productRewards, acc);
		rewardsProduct = StateValues.createFromDoubleArray(res.soln, product.getProductModel());

		// Output vector over product, if required
		if (getExportProductVector()) {
			mainLog.println("\nExporting product solution vector matrix to file \"" + getExportProductVectorFilename() + "\"...");
			PrismFileLog out = new PrismFileLog(getExportProductVectorFilename());
			rewardsProduct.print(out, false, false, false, false);
			out.close();
		}

		// Mapping rewards in the original model
		rewards = product.projectToOriginalModel(rewardsProduct);
		rewardsProduct.clear();

		return rewards;
	}

	public ModelCheckerResult computeInstantaneousRewards(DTMC dtmc, MCRewards mcRewards, double t) throws PrismException
	{
		ModelCheckerResult res;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		long timer;
		int right = (int) t;

		// Store num states
		n = dtmc.getNumStates();

		// Start backwards transient computation
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting backwards instantaneous rewards computation...");

		// Create solution vector(s)
		soln = new double[n];
		soln2 = new double[n];

		// Initialise solution vectors.
		for (i = 0; i < n; i++)
			soln[i] = mcRewards.getStateReward(i);

		// Start iterations
		for (iters = 0; iters < right; iters++) {
			// Matrix-vector multiply
			dtmc.mvMult(soln, soln2, null);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished backwards transient computation
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Backwards transient instantaneous rewards computation");
		mainLog.println(" took " + iters + " iters and " + timer / 1000.0 + " seconds.");

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.lastSoln = soln2;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		res.timePre = 0.0;
		return res;
	}

	public ModelCheckerResult computeCumulativeRewards(DTMC dtmc, MCRewards mcRewards, double t) throws PrismException
	{
		ModelCheckerResult res;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		long timer;
		int right = (int) t;

		// Store num states
		n = dtmc.getNumStates();

		// Start backwards transient computation
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting backwards cumulative rewards computation...");

		// Create solution vector(s)
		soln = new double[n];
		soln2 = new double[n];

		// Start iterations
		for (iters = 0; iters < right; iters++) {
			// Matrix-vector multiply plus adding rewards
			dtmc.mvMult(soln, soln2, null);
			for (i = 0; i < n; i++) {
				soln2[i] += mcRewards.getStateReward(i);
			}
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished backwards transient computation
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Backwards cumulative rewards computation");
		mainLog.println(" took " + iters + " iters and " + timer / 1000.0 + " seconds.");

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.lastSoln = soln2;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		res.timePre = 0.0;
		return res;
	}

	public ModelCheckerResult computeTotalRewards(DTMC dtmc, MCRewards mcRewards) throws PrismException
	{
		ModelCheckerResult res;
		int n, numBSCCs;
		long timer;

		// Switch to a supported method, if necessary
		if (!(linEqMethod == LinEqMethod.POWER)) {
			linEqMethod = LinEqMethod.POWER;
			mainLog.printWarning("Switching to linear equation solution method \"" + linEqMethod.fullName() + "\"");
		}

		// Store num states
		n = dtmc.getNumStates();

		// Start total rewards computation
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting total reward computation...");

		// Compute bottom strongly connected components (BSCCs)
		SCCComputer sccComputer = SCCComputer.createSCCComputer(this, dtmc);
		sccComputer.computeBSCCs();
		List<NatBitSet> bsccs = sccComputer.getBSCCs();
		numBSCCs = bsccs.size();

		// Find BSCCs with non-zero reward
		NatBitSet bsccsNonZero = NatBitSets.boundedSet(n);
		for (int b = 0; b < numBSCCs; b++) {
			NatBitSet bscc = bsccs.get(b);
			IntIterator iterator = bscc.iterator();
			while (iterator.hasNext()) {
				int i = iterator.nextInt();
				if (mcRewards.getStateReward(i) > 0) {
					bsccsNonZero.or(bscc);
					break;
				}
			}
		}
		mainLog.print("States in non-zero reward BSCCs: " + bsccsNonZero.size());

		// Find states with infinite reward (those reach a non-zero reward BSCC with prob > 0)
		NatBitSet inf = prob0(dtmc, null, bsccsNonZero).complement();
		int numInf = inf.size();
		mainLog.println(", inf=" + numInf + ", maybe=" + (n - numInf));

		// Compute rewards
		// (do this using the functions for "reward reachability" properties but with no targets)
		switch (linEqMethod) {
		case POWER:
			res = computeReachRewardsValIter(dtmc, mcRewards, NatBitSets.boundedSet(dtmc.getNumStates()), inf, null, null);
			break;
		default:
			throw new PrismException("Unknown linear equation solution method " + linEqMethod.fullName());
		}

		// Finished total reward computation
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Total reward computation");
		mainLog.println(" took " + timer / 1000.0 + " seconds.");

		// Return results
		return res;
	}

	// Steady-state/transient probability computation

	/**
	 * Compute steady-state probability distribution (forwards).
	 * Start from initial state (or uniform distribution over multiple initial states).
	 */
	public StateValues doSteadyState(DTMC dtmc) throws PrismException
	{
		return doSteadyState(dtmc, (StateValues) null);
	}

	/**
	 * Compute steady-state probability distribution (forwards).
	 * Optionally, use the passed in file initDistFile to give the initial probability distribution (time 0).
	 * If null, start from initial state (or uniform distribution over multiple initial states).
	 */
	public StateValues doSteadyState(DTMC dtmc, File initDistFile) throws PrismException
	{
		StateValues initDist = readDistributionFromFile(initDistFile, dtmc);
		return doSteadyState(dtmc, initDist);
	}

	/**
	 * Compute steady-state probability distribution (forwards).
	 * Optionally, use the passed in vector initDist as the initial probability distribution (time 0).
	 * If null, start from initial state (or uniform distribution over multiple initial states).
	 * For reasons of efficiency, when a vector is passed in, it will be trampled over,
	 * so if you wanted it, take a copy.
	 *
	 * @param dtmc     The DTMC
	 * @param initDist Initial distribution (will be overwritten)
	 */
	public StateValues doSteadyState(DTMC dtmc, StateValues initDist) throws PrismException
	{
		StateValues initDistNew = (initDist == null) ? buildInitialDistribution(dtmc) : initDist;
		ModelCheckerResult res = computeSteadyStateProbs(dtmc, initDistNew.getDoubleArray());
		return StateValues.createFromDoubleArray(res.soln, dtmc);
	}

	/**
	 * Compute transient probability distribution (forwards).
	 * Optionally, use the passed in vector initDist as the initial probability distribution (time step 0).
	 * If null, start from initial state (or uniform distribution over multiple initial states).
	 * For reasons of efficiency, when a vector is passed in, it will be trampled over,
	 * so if you wanted it, take a copy.
	 *
	 * @param dtmc     The DTMC
	 * @param k        Time step
	 * @param initDist Initial distribution (will be overwritten)
	 */
	public StateValues doTransient(DTMC dtmc, int k, double initDist[]) throws PrismException
	{
		throw new PrismNotSupportedException("Not implemented yet");
	}

	// Numerical computation functions

	/**
	 * Compute next=state probabilities.
	 * i.e. compute the probability of being in a state in {@code target} in the next step.
	 *
	 * @param dtmc   The DTMC
	 * @param target Target states
	 */
	public ModelCheckerResult computeNextProbs(DTMC dtmc, NatBitSet target) throws PrismException
	{
		ModelCheckerResult res;
		int n;
		double soln[], soln2[];
		long timer;

		timer = System.currentTimeMillis();

		// Store num states
		n = dtmc.getNumStates();

		// Create/initialise solution vector(s)
		soln = Utils.intSetToDoubleArray(target, n);
		soln2 = new double[n];

		// Next-step probabilities
		dtmc.mvMult(soln, soln2, null);

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln2;
		res.numIters = 1;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Given a value vector x, compute the probability:
	 * v(s) = Sum_s' P(s,s')*x(s')   for s labeled with a,
	 * v(s) = 0                      for s not labeled with a.
	 *
	 * @param dtmc the DTMC model
	 * @param a    the set of states labeled with a
	 * @param x    the value vector
	 */
	protected double[] computeRestrictedNext(DTMC dtmc, NatBitSet a, double[] x)
	{
		double[] soln;
		int n;

		// Store num states
		n = dtmc.getNumStates();

		// initialized to 0.0
		soln = new double[n];

		// Next-step probabilities multiplication
		// restricted to a states
		dtmc.mvMult(x, soln, a);

		return soln;
	}

	/**
	 * Compute reachability probabilities.
	 * i.e. compute the probability of reaching a state in {@code target}.
	 *
	 * @param dtmc   The DTMC
	 * @param target Target states
	 */
	public ModelCheckerResult computeReachProbs(DTMC dtmc, NatBitSet target) throws PrismException
	{
		return computeReachProbs(dtmc, null, target, null, null);
	}

	/**
	 * Compute until probabilities.
	 * i.e. compute the probability of reaching a state in {@code target},
	 * while remaining in those in {@code remain}.
	 *
	 * @param dtmc   The DTMC
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 */
	public ModelCheckerResult computeUntilProbs(DTMC dtmc, NatBitSet remain, NatBitSet target) throws PrismException
	{
		return computeReachProbs(dtmc, remain, target, null, null);
	}

	/**
	 * Compute reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in {@code remain}.
	 *
	 * @param dtmc   The DTMC
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param init   Optionally, an initial solution vector (may be overwritten)
	 * @param known  Optionally, a set of states for which the exact answer is known
	 *               Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	public ModelCheckerResult computeReachProbs(DTMC dtmc, NatBitSet remain, NatBitSet target, double init[], NatBitSet known) throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet no, yes;
		int n, numYes, numNo;
		long timer, timerProb0, timerProb1;
		PredecessorRelation pre = null;
		// Local copy of setting
		LinEqMethod linEqMethod = this.linEqMethod;
		DTMCSolnMethod dtmcSolnMethod = this.dtmcSolnMethod;

		if (dtmcSolnMethod == DTMCSolnMethod.LINEAR_EQUATIONS) {
			// Switch to a supported method, if necessary
			if (!EnumSet.of(LinEqMethod.GAUSS_SEIDEL, LinEqMethod.POWER, LinEqMethod.DECOMPOSITION).contains(linEqMethod)) {
				linEqMethod = LinEqMethod.GAUSS_SEIDEL;
				mainLog.printWarning("Switching to linear equation solution method \"" + linEqMethod.fullName() + "\"");
			}
		}

		// Start probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting probabilistic reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		dtmc.checkForDeadlocks(target);

		// Store num states
		n = dtmc.getNumStates();

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

		// If required, export info about target states
		if (getExportTarget()) {
			NatBitSet bsInit = NatBitSets.boundedSet(n);
			for (int i = 0; i < n; i++) {
				bsInit.set(i, dtmc.isInitialState(i));
			}
			List<NatBitSet> labels = Arrays.asList(bsInit, target);
			List<String> labelNames = Arrays.asList("init", "target");
			mainLog.println("\nExporting target states info to file \"" + getExportTargetFilename() + "\"...");
			exportLabels(dtmc, labels, labelNames, Prism.EXPORT_PLAIN, new PrismFileLog(getExportTargetFilename()));
		}

		if (precomp && (prob0 || prob1) && preRel) {
			pre = dtmc.getPredecessorRelation(this, true);
		}

		// Precomputation
		timerProb0 = System.currentTimeMillis();
		if (precomp && prob0) {
			if (preRel) {
				no = prob0(dtmc, remain, target, pre);
			} else {
				no = prob0(dtmc, remain, target);
			}
		} else {
			no = NatBitSets.emptySet();
		}
		timerProb0 = System.currentTimeMillis() - timerProb0;
		timerProb1 = System.currentTimeMillis();
		if (precomp && prob1) {
			if (preRel) {
				yes = prob1(dtmc, remain, target, pre);
			} else {
				yes = prob1(dtmc, remain, target);
			}
		} else {
			yes = target.clone();
		}
		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numYes = yes.size();
		numNo = no.size();
		mainLog.println(String.format("target=%d, yes=%d, no=%d, maybe=%d", target.size(), numYes, numNo, n - (numYes + numNo)));

		// Compute probabilities
		if (dtmcSolnMethod == DTMCSolnMethod.LINEAR_EQUATIONS) {
			switch (linEqMethod) {
			case DECOMPOSITION:
				res = computeReachProbsLinearEquations(dtmc, no, yes, init == null ? null : s -> init[s], known);
				break;
			case POWER:
				res = computeReachProbsValIter(dtmc, no, yes, init, known);
				break;
			case GAUSS_SEIDEL:
				res = computeReachProbsGaussSeidel(dtmc, no, yes, init, known);
				break;
			default:
				throw new PrismException("Unknown linear equation solution method " + linEqMethod.fullName());
			}
		} else if (dtmcSolnMethod == DTMCSolnMethod.INTERVAL_ITERATION) {
			ProcessingResult<IntervalResultDense> result = computeReachProbsIntervalValIter(dtmc, no, yes, init, init, known, termCritParam);
			IntervalResultDense intervalResult = result.getResult();
			res = new ModelCheckerResult();
			res.numIters = result.getNumIters();
			res.soln = new double[n];
			for (int i = 0; i < n; i++) {
				res.soln[i] = intervalResult.getAverage(i);
			}
		} else {
			throw new PrismException("Unknown dtmc solution method " + dtmcSolnMethod.fullName());
		}

		// Finished probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Probabilistic reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		res.timeProb0 = timerProb0 / 1000.0;
		res.timePre = (timerProb0 + timerProb1) / 1000.0;

		return res;
	}

	/**
	 * Prob0 precomputation algorithm (using predecessor relation),
	 * i.e. determine the states of a DTMC which, with probability 0,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 *
	 * @param dtmc   The DTMC
	 * @param remain Remain in these states (optional: {@code null} means "all states")
	 * @param target Target states
	 * @param pre    The predecessor relation
	 */
	public NatBitSet prob0(DTMC dtmc, NatBitSet remain, NatBitSet target, PredecessorRelation pre)
	{
		NatBitSet canReachTarget, result;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		mainLog.println("Starting Prob0...");

		// Special case: no target states
		if (target.isEmpty()) {
			return NatBitSets.boundedFilledSet(dtmc.getNumStates());
		}

		// calculate all states that can reach 'target'
		// while remaining in 'remain' in the underlying graph,
		// where all the 'target' states are made absorbing
		canReachTarget = PredecessorRelations.calculatePreStar(pre, remain, target, target);

		// prob0 = complement of 'canReachTarget'
		result = NatBitSets.ensureBounded(canReachTarget, dtmc.getNumStates()).complement();

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Prob0");
		mainLog.println(" took " + timer / 1000.0 + " seconds.");

		return result;
	}

	/**
	 * Prob0 precomputation algorithm (using a fixed-point computation),
	 * i.e. determine the states of a DTMC which, with probability 0,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 *
	 * @param dtmc   The DTMC
	 * @param remain Remain in these states (optional: {@code null} means "all")
	 * @param target Target states
	 */
	public BoundedNatBitSet prob0(DTMC dtmc, NatBitSet remain, NatBitSet target)
	{
		int n, iters;
		BoundedNatBitSet u, soln, unknown;
		boolean u_done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		mainLog.println("Starting Prob0...");

		// Special case: no target states
		if (target.isEmpty()) {
			return NatBitSets.boundedFilledSet(dtmc.getNumStates());
		}

		// Initialise vectors
		n = dtmc.getNumStates();
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
			dtmc.prob0step(unknown, u, soln);
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
		mainLog.print("Prob0");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		return u;
	}

	/**
	 * Prob1 precomputation algorithm (using predecessor relation),
	 * i.e. determine the states of a DTMC which, with probability 1,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 *
	 * @param dtmc   The DTMC
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param pre    The predecessor relation of the DTMC
	 */
	public NatBitSet prob1(DTMC dtmc, NatBitSet remain, NatBitSet target, PredecessorRelation pre)
	{
		// Implements the constrained reachability algorithm from
		// Baier, Katoen: Principles of Model Checking (Corollary 10.31 Qualitative Constrained Reachability)
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		mainLog.println("Starting Prob1...", PrismLog.VL_ALL);

		// Special case: no 'target' states
		if (target.isEmpty()) {
			return NatBitSets.emptySet();
		}

		// mark all states in 'target' and all states not in 'remain' as absorbing
		NatBitSet absorbing;
		if (remain == null) {
			// for remain == null, remain consists of all states, thus, absorbing = the empty set is already the complement of remain
			absorbing = NatBitSets.modifiableCopyOf(target);
		} else {
			// complement remain and union with target
			absorbing = NatBitSets.ensureBounded(remain, dtmc.getNumStates()).complement();
			absorbing.or(target);
		}

		// M' = DTMC where all 'absorbing' states are considered to be absorbing

		// the set of states that satisfy E [ F target ] in M'
		// Pre*(target)
		NatBitSet canReachTarget = PredecessorRelations.calculatePreStar(pre, null, target, absorbing);
		assert remain == null || Sets.union(target, remain).containsAll(canReachTarget);

		// complement canReachTarget
		// S\Pre*(target)
		NatBitSet canNotReachTarget = NatBitSets.asBounded(canReachTarget, dtmc.getNumStates()).complement();

		// the set of states that can reach a canNotReachTarget state in M'
		// Pre*(S\Pre*(target))
		NatBitSet probTargetNot1 = PredecessorRelations.calculatePreStar(pre, null, canNotReachTarget, absorbing);

		// complement probTargetNot1
		// S\Pre*(S\Pre*(target))
		NatBitSet result = NatBitSets.asBounded(probTargetNot1, dtmc.getNumStates()).complement();

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Prob1", PrismLog.VL_ALL);
		mainLog.println(" took " + timer / 1000.0 + " seconds.", PrismLog.VL_ALL);

		return result;
	}

	/**
	 * Prob1 precomputation algorithm (using a fixed-point computation)
	 * i.e. determine the states of a DTMC which, with probability 1,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 *
	 * @param dtmc   The DTMC
	 * @param remain Remain in these states (optional: {@code null} means "all")
	 * @param target Target states
	 */
	public NatBitSet prob1(DTMC dtmc, NatBitSet remain, NatBitSet target)
	{
		int n, iters;
		NatBitSet u, v, soln, unknown;
		boolean u_done, v_done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		mainLog.println("Starting Prob1...");

		// Special case: no target states
		if (target.isEmpty()) {
			return NatBitSets.setWithExpectedSize(dtmc.getNumStates());
		}

		// Initialise vectors
		n = dtmc.getNumStates();
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
				dtmc.prob1step(unknown, u, v, soln);
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
		mainLog.print("Prob1");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		return u;
	}

	/**
	 * Interval iteration for DTMCs inspired by
	 * Serge Haddad, Benjamin Monmege, Interval iteration algorithm for MDPs and IMDPs, Theoretical Computer Science,
	 * Available online 2 January 2017, ISSN 0304-3975, http://dx.doi.org/10.1016/j.tcs.2016.12.003.
	 * <p>
	 * Maintains upper (initialized to 1) and lower bounds (intialized to 0) for all states and performs
	 * value iteration on them. The update is performed both over the upper as well as the lower values.
	 * The difference between the upper and lower values gives the error bound.
	 */
	public ProcessingResult<IntervalResultDense> computeReachProbsIntervalValIter(DTMC dtmc, NatBitSet no, NatBitSet yes, double initLower[],
			double initUpper[],
			NatBitSet known, double precision) throws PrismException
	{
		Time time = new Time();
		mainLog.println("Starting reachability interval value iteration...");

		int n = dtmc.getNumStates();
		double[] lower = new double[n];
		double[] upper = new double[n];
		double[] lower2 = new double[n];
		double[] upper2 = new double[n];
		Arrays.fill(upper, 1d);
		Arrays.fill(upper2, 1d);

		IntIterator noIter = no.iterator();
		while (noIter.hasNext()) {
			int state = noIter.nextInt();
			upper[state] = 0d;
			upper2[state] = 0d;
		}
		IntIterator yesIter = yes.iterator();
		while (yesIter.hasNext()) {
			int state = yesIter.nextInt();
			lower[state] = 1d;
			lower2[state] = 1d;
		}
		if (known != null) {
			IntIterator knownIter = known.iterator();
			while (knownIter.hasNext()) {
				int knownState = knownIter.nextInt();
				assert initLower[knownState] <= initUpper[knownState];
				lower[knownState] = initLower[knownState];
				upper[knownState] = initUpper[knownState];
			}
		}

		NatBitSet unknown = NatBitSets.boundedFilledSet(n);
		unknown.andNot(yes);
		unknown.andNot(no);
		if (known != null) {
			unknown.andNot(known);
		}

		// Start iterations
		int iterations = 0;
		boolean done = false;
		while (!done && iterations < maxIters) {
			iterations++;

			dtmc.mvMult(lower, lower2, unknown);
			double[] swapLower = lower;
			lower = lower2;
			lower2 = swapLower;
			dtmc.mvMult(upper, upper2, unknown);
			double[] swapUpper = upper;
			upper = upper2;
			upper2 = swapUpper;

			done = true;
			for (int state = 0; state < upper.length; state++) {
				if (!PrismUtils.doublesAreCloseAbs(lower[state], upper[state], precision)) {
					done = false;
					// unknown.clear(state);
					break;
				}
			}
		}

		// Finished value iteration
		double elapsedSeconds = time.elapsedSeconds();
		mainLog.println(String.format("Value iteration took %d iterations and %f seconds", iterations, elapsedSeconds));

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iterations + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Return results
		return new ProcessingResult<>(new IntervalResultDense(upper, lower), precision, iterations, elapsedSeconds);
	}

	/**
	 * Compute reachability probabilities using value iteration.
	 *
	 * @param dtmc  The DTMC
	 * @param no    Probability 0 states
	 * @param yes   Probability 1 states
	 * @param init  Optionally, an initial solution vector (will be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 *              Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachProbsValIter(DTMC dtmc, NatBitSet no, NatBitSet yes, double init[], NatBitSet known) throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[], initVal;
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting value iteration...");

		// Store num states
		n = dtmc.getNumStates();

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

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			iters++;
			// Matrix-vector multiply
			dtmc.mvMult(soln, soln2, unknown);
			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Value iteration");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

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

	public ModelCheckerResult computeReachProbsLinearEquations(DTMC dtmc, NatBitSet no, NatBitSet yes, IntToDoubleFunction init, NatBitSet known)
			throws PrismException
	{
		checkArgument(!yes.intersects(no));

		Time time = new Time();

		if (known != null && (known.intersects(yes) || known.intersects(no))) {
			known = known.clone();
			known.andNot(yes);
			known.andNot(no);
		}

		int numDtmcStates = dtmc.getNumStates();
		int numNoStates = no == null ? 0 : no.size();
		int numYesStates = yes.size();
		int numKnownStates = known == null ? 0 : known.size();

		int numUnknownStates = numDtmcStates - numNoStates - numYesStates - numKnownStates;
		NatBitSet unknownStates = NatBitSets.boundedFilledSet(numUnknownStates);
		unknownStates.andNot(yes);
		if (no != null) {
			unknownStates.andNot(no);
		}
		if (known != null) {
			unknownStates.and(known);
		}

		SparseStore<Double> systemLHS = SparseStore2.makePrimitive(numUnknownStates, numUnknownStates);
		PrimitiveDenseStore systemRHS = PrimitiveDenseStore.FACTORY.makeZero(numUnknownStates, 1);

		// Solve reach(s) = sum_{s' in unknown} p(s, s') * reach(s') + sum_{s' \in target} p(s, s')

		IntUnaryOperator stateIndexMap = FastUtils.elementToIndexMap(unknownStates);
		IntIterator stateIterator = unknownStates.iterator();
		while (stateIterator.hasNext()) {
			int unknownState = stateIterator.nextInt();
			int unknownStateIndex = stateIndexMap.applyAsInt(unknownState);

			Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(unknownState);
			while (transitionIterator.hasNext()) {
				Int2DoubleMap.Entry transition = transitionIterator.next();
				int destinationState = transition.getIntKey();
				double prob = transition.getDoubleValue();

				if (yes.contains(destinationState)) {
					systemRHS.add(unknownStateIndex, prob);
				} else if (known != null && known.contains(destinationState)) {
					systemRHS.add(unknownStateIndex, prob * init.applyAsDouble(destinationState));
				} else if (no == null || !no.contains(destinationState)) {
					int destinationStateIndex = stateIndexMap.applyAsInt(destinationState);
					// P
					systemLHS.set(unknownStateIndex, destinationStateIndex, prob);
				}
			}
			// -I
			systemLHS.add(unknownStateIndex, unknownStateIndex, -1.0d);
		}

		MatrixStore<Double> solution = solveLinearEquationSystem(systemLHS, systemRHS, false);

		double[] soln = new double[numDtmcStates];

		yes.forEach((IntConsumer) s -> soln[s] = 1.0d);
		if (no != null) {
			no.forEach((IntConsumer) s -> soln[s] = 0.0d);
		}
		if (known != null) {
			known.forEach((IntConsumer) s -> soln[s] = init.applyAsDouble(s));
		}
		stateIterator = unknownStates.iterator();
		while (stateIterator.hasNext()) {
			int state = stateIterator.nextInt();
			int stateIndex = stateIndexMap.applyAsInt(state);
			soln[state] = solution.doubleValue(stateIndex);
		}

		ModelCheckerResult res = new ModelCheckerResult();
		res.soln = soln;
		res.timeTaken = time.elapsedSeconds();
		return res;
	}

	/**
	 * Compute reachability probabilities using Gauss-Seidel.
	 *
	 * @param dtmc  The DTMC
	 * @param no    Probability 0 states
	 * @param yes   Probability 1 states
	 * @param init  Optionally, an initial solution vector (will be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 *              Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	public ModelCheckerResult computeReachProbsGaussSeidel(DTMC dtmc, NatBitSet no, NatBitSet yes, double init[], NatBitSet known) throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], initVal, maxDiff;
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting Gauss-Seidel...");

		// Store num states
		n = dtmc.getNumStates();

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
			// Matrix-vector multiply
			maxDiff = dtmc.mvMultGS(soln, unknown, termCrit == TermCrit.ABSOLUTE);
			// Check termination
			done = maxDiff < termCritParam;
		}

		// Finished Gauss-Seidel
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Gauss-Seidel");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

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
	 * Compute bounded reachability probabilities.
	 * i.e. compute the probability of reaching a state in {@code target} within k steps.
	 *
	 * @param dtmc   The DTMC
	 * @param target Target states
	 * @param k      Bound
	 */
	public ModelCheckerResult computeBoundedReachProbs(DTMC dtmc, NatBitSet target, int k) throws PrismException
	{
		return computeBoundedReachProbs(dtmc, null, target, k, null, null);
	}

	/**
	 * Compute bounded until probabilities.
	 * i.e. compute the probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 *
	 * @param dtmc   The DTMC
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k      Bound
	 */
	public ModelCheckerResult computeBoundedUntilProbs(DTMC dtmc, NatBitSet remain, NatBitSet target, int k) throws PrismException
	{
		return computeBoundedReachProbs(dtmc, remain, target, k, null, null);
	}

	/**
	 * Compute bounded reachability/until probabilities.
	 * i.e. compute the probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 *
	 * @param dtmc    The DTMC
	 * @param remain  Remain in these states (optional: null means "all")
	 * @param target  Target states
	 * @param k       Bound
	 * @param init    Initial solution vector - pass null for default
	 * @param results Optional array of size b+1 to store (init state) results for each step (null if unused)
	 */
	public ModelCheckerResult computeBoundedReachProbs(DTMC dtmc, NatBitSet remain, NatBitSet target, int k, double init[], double results[])
			throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		long timer;

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting bounded probabilistic reachability...");

		// Store num states
		n = dtmc.getNumStates();

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
			// TODO: whether this is min or max should be specified somehow
			results[0] = Utils.minMaxOverArraySubset(soln2, dtmc.getInitialStates(), true);
		}

		// Determine set of states actually need to perform computation for
		unknown = NatBitSets.boundedFilledSet(n);
		unknown.andNot(target);
		if (remain != null)
			unknown.and(remain);

		// Start iterations
		iters = 0;
		while (iters < k) {

			iters++;
			// Matrix-vector multiply
			dtmc.mvMult(soln, soln2, unknown);
			// Store intermediate results if required
			// (compute min/max value over initial states for this step)
			if (results != null) {
				// TODO: whether this is min or max should be specified somehow
				results[iters] = Utils.minMaxOverArraySubset(soln2, dtmc.getInitialStates(), true);
			}
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Bounded probabilistic reachability");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

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
	 *
	 * @param dtmc      The DTMC
	 * @param mcRewards The rewards
	 * @param target    Target states
	 */
	public ModelCheckerResult computeReachRewards(DTMC dtmc, MCRewards mcRewards, NatBitSet target) throws PrismException
	{
		return computeReachRewards(dtmc, mcRewards, target, null, null);
	}

	/**
	 * Compute expected reachability rewards.
	 *
	 * @param dtmc      The DTMC
	 * @param mcRewards The rewards
	 * @param target    Target states
	 * @param init      Optionally, an initial solution vector (may be overwritten)
	 * @param known     Optionally, a set of states for which the exact answer is known
	 *                  Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	public ModelCheckerResult computeReachRewards(DTMC dtmc, MCRewards mcRewards, NatBitSet target, double init[], NatBitSet known) throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet inf;
		int n, numTarget, numInf;
		long timer, timerProb1;
		// Local copy of setting
		LinEqMethod linEqMethod = this.linEqMethod;

		// Switch to a supported method, if necessary
		if (!(linEqMethod == LinEqMethod.POWER)) {
			linEqMethod = LinEqMethod.POWER;
			mainLog.printWarning("Switching to linear equation solution method \"" + linEqMethod.fullName() + "\"");
		}

		// Start expected reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting expected reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		dtmc.checkForDeadlocks(target);

		// Store num states
		n = dtmc.getNumStates();

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
		inf = prob1(dtmc, null, target);
		inf.flip(0, n);
		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numTarget = target.size();
		numInf = inf.size();
		mainLog.println("target=" + numTarget + ", inf=" + numInf + ", rest=" + (n - (numTarget + numInf)));

		// Compute rewards
		switch (linEqMethod) {
		case POWER:
			res = computeReachRewardsValIter(dtmc, mcRewards, target, inf, init, known);
			break;
		default:
			throw new PrismException("Unknown linear equation solution method " + linEqMethod.fullName());
		}

		// Finished expected reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Expected reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		res.timePre = timerProb1 / 1000.0;

		return res;
	}

	/**
	 * Compute expected reachability rewards using value iteration.
	 *
	 * @param dtmc      The DTMC
	 * @param mcRewards The rewards
	 * @param target    Target states
	 * @param inf       States for which reward is infinite
	 * @param init      Optionally, an initial solution vector (will be overwritten)
	 * @param known     Optionally, a set of states for which the exact answer is known
	 *                  Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	public ModelCheckerResult computeReachRewardsValIter(DTMC dtmc, MCRewards mcRewards, NatBitSet target, NatBitSet inf, double init[], NatBitSet known)
			throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting value iteration...");

		// Store num states
		n = dtmc.getNumStates();

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
		if (known != null) {
			unknown.andNot(known);
		}

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			//mainLog.println(soln);
			iters++;
			// Matrix-vector multiply
			dtmc.mvMult(soln, soln2, unknown);
			IntIterator iterator = unknown.iterator();
			while (iterator.hasNext()) {
				int state = iterator.nextInt();
				soln2[state] += mcRewards.getStateReward(state);
			}

			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Value iteration");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

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
	 * Compute (forwards) steady-state probabilities
	 * i.e. compute the long-run probability of being in each state,
	 * assuming the initial distribution {@code initDist}.
	 * For space efficiency, the initial distribution vector will be modified and values over-written,
	 * so if you wanted it, take a copy.
	 *
	 * @param dtmc     The DTMC
	 * @param initDist Initial distribution (will be overwritten)
	 */
	public ModelCheckerResult computeSteadyStateProbs(DTMC dtmc, double initDist[]) throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet startNot, bscc;
		double probBSCCs[], solnProbs[], reachProbs[];
		int n, numBSCCs, allInOneBSCC;
		long timer;

		timer = System.currentTimeMillis();

		// Store num states
		n = dtmc.getNumStates();
		// Create results vector
		solnProbs = new double[n];

		// Compute bottom strongly connected components (BSCCs)
		SCCComputer sccComputer = SCCComputer.createSCCComputer(this, dtmc);
		sccComputer.computeBSCCs();
		List<NatBitSet> bsccs = sccComputer.getBSCCs();
		NatBitSet notInBSCCs = sccComputer.getNotInBSCCs();
		numBSCCs = bsccs.size();

		// See which states in the initial distribution do *not* have non-zero prob
		startNot = NatBitSets.boundedSet(n);
		for (int i = 0; i < n; i++) {
			if (initDist[i] == 0)
				startNot.set(i);
		}
		// Determine whether initial states are all in a single BSCC
		allInOneBSCC = -1;
		for (int b = 0; b < numBSCCs; b++) {
			if (!bsccs.get(b).intersects(startNot)) {
				allInOneBSCC = b;
				break;
			}
		}

		// If all initial states are in a single BSCC, it's easy...
		// Just compute steady-state probabilities for the BSCC
		if (allInOneBSCC != -1) {
			mainLog.println("\nInitial states all in one BSCC (so no reachability probabilities computed)");
			bscc = bsccs.get(allInOneBSCC);
			computeSteadyStateProbsForBSCC(dtmc, bscc, solnProbs);
		}

		// Otherwise, have to consider all the BSCCs
		else {

			// Compute probability of reaching each BSCC from initial distribution
			probBSCCs = new double[numBSCCs];
			for (int b = 0; b < numBSCCs; b++) {
				mainLog.println("\nComputing probability of reaching BSCC " + (b + 1));
				bscc = bsccs.get(b);
				// Compute probabilities
				reachProbs = computeUntilProbs(dtmc, notInBSCCs, bscc).soln;
				// Compute probability of reaching BSCC, which is dot product of
				// vectors for initial distribution and probabilities of reaching it
				probBSCCs[b] = 0.0;
				for (int i = 0; i < n; i++) {
					probBSCCs[b] += initDist[i] * reachProbs[i];
				}
				mainLog.print("\nProbability of reaching BSCC " + (b + 1) + ": " + probBSCCs[b] + "\n");
			}

			// Compute steady-state probabilities for each BSCC
			for (int b = 0; b < numBSCCs; b++) {
				mainLog.println("\nComputing steady-state probabilities for BSCC " + (b + 1));
				bscc = bsccs.get(b);
				// Compute steady-state probabilities for the BSCC
				computeSteadyStateProbsForBSCC(dtmc, bscc, solnProbs);
				// Multiply by BSCC reach prob
				double probBSCC = probBSCCs[b];
				bscc.forEach((IntConsumer) i -> solnProbs[i] *= probBSCC);
			}
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = solnProbs;
		timer = System.currentTimeMillis() - timer;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Perform (backwards) steady-state probabilities, as required for (e.g. CSL) model checking.
	 * Compute, for each initial state s, the sum over all states s'
	 * of the steady-state probability of being in s'
	 * multiplied by the corresponding probability in the vector {@code multProbs}.
	 * If {@code multProbs} is null, it is assumed to be all 1s.
	 *
	 * @param dtmc      The DTMC
	 * @param multProbs Multiplication vector (optional: null means all 1s)
	 */
	public ModelCheckerResult computeSteadyStateBackwardsProbs(DTMC dtmc, double multProbs[]) throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet bscc;
		double probBSCCs[], ssProbs[], reachProbs[], soln[];
		int n, numBSCCs;
		long timer;

		timer = System.currentTimeMillis();

		// Store num states
		n = dtmc.getNumStates();

		// Compute bottom strongly connected components (BSCCs)
		SCCComputer sccComputer = SCCComputer.createSCCComputer(this, dtmc);
		sccComputer.computeBSCCs();
		List<NatBitSet> bsccs = sccComputer.getBSCCs();
		NatBitSet notInBSCCs = sccComputer.getNotInBSCCs();
		numBSCCs = bsccs.size();

		// Compute steady-state probability for each BSCC...
		probBSCCs = new double[numBSCCs];
		ssProbs = new double[n];
		for (int b = 0; b < numBSCCs; b++) {
			mainLog.println("\nComputing steady state probabilities for BSCC " + (b + 1));
			bscc = bsccs.get(b);
			// Compute steady-state probabilities for the BSCC
			computeSteadyStateProbsForBSCC(dtmc, bscc, ssProbs);
			// Compute weighted sum of probabilities with multProbs
			probBSCCs[b] = 0.0;
			IntIterator iterator = bscc.iterator();
			if (multProbs == null) {
				while (iterator.hasNext()) {
					probBSCCs[b] += ssProbs[iterator.nextInt()];
				}
			} else {
				while (iterator.hasNext()) {
					int i = iterator.nextInt();
					probBSCCs[b] += multProbs[i] * ssProbs[i];
				}
			}
			mainLog.print("\nValue for BSCC " + (b + 1) + ": " + probBSCCs[b] + "\n");
		}

		// Create/initialise prob vector
		soln = new double[n];
		for (int i = 0; i < n; i++) {
			soln[i] = 0.0;
		}

		// If every state is in a BSCC, it's much easier...
		if (notInBSCCs.isEmpty()) {
			mainLog.println("\nAll states are in BSCCs (so no reachability probabilities computed)");
			for (int b = 0; b < numBSCCs; b++) {
				bscc = bsccs.get(b);
				double probBSCC = probBSCCs[b];
				bscc.forEach((IntConsumer) i -> soln[i] += probBSCC);
			}
		}

		// Otherwise we have to do more work...
		else {
			// Compute probabilities of reaching each BSCC...
			for (int b = 0; b < numBSCCs; b++) {
				// Skip BSCCs with zero probability
				if (probBSCCs[b] == 0.0) {
					continue;
				}
				mainLog.println("\nComputing probabilities of reaching BSCC " + (b + 1));
				bscc = bsccs.get(b);
				// Compute probabilities
				reachProbs = computeUntilProbs(dtmc, notInBSCCs, bscc).soln;
				// Multiply by value for BSCC, add to total
				for (int i = 0; i < n; i++) {
					soln[i] += reachProbs[i] * probBSCCs[b];
				}
			}
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		timer = System.currentTimeMillis() - timer;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute steady-state probabilities for a BSCC
	 * i.e. compute the long-run probability of being in each state of the BSCC.
	 * No initial distribution is specified since it does not affect the result.
	 * The result will be stored in the relevant portion of a full vector,
	 * whose size equals the number of states in the DTMC.
	 * Optionally, pass in an existing vector to be used for this purpose.
	 *
	 * @param dtmc   The DTMC
	 * @param bscc   The BSCC to be analysed
	 * @param result Storage for result (ignored if null)
	 */
	public ModelCheckerResult computeSteadyStateProbsForBSCC(DTMC dtmc, NatBitSet bscc, double result[]) throws PrismException
	{
		ModelCheckerResult res;
		int n, iters;
		double soln[], soln2[], tmpsoln[];
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting value iteration...");

		// Store num states
		n = dtmc.getNumStates();

		// Create solution vector(s)
		// Use the passed in vector, if present
		soln = result == null ? new double[n] : result;
		soln2 = new double[n];

		// Initialise solution vectors. Equiprobable for BSCC states.
		double equiprob = 1.0 / bscc.size();
		IntIterator iterator = bscc.iterator();
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			soln[i] = equiprob;
			soln2[i] = equiprob;
		}

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			iters++;
			// Matrix-vector multiply
			dtmc.vmMult(soln, soln2);
			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Value iteration");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

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

	private MatrixStore<Double> solveLinearEquationSystem(MatrixStore<Double> systemLHS, MatrixStore<Double> systemRHS, boolean useSvd)
			throws PrismException
	{
		SolverTask<Double> solver;
		if (useSvd) {
			solver = SingularValue.make(systemLHS);
		} else {
			solver = LU.make(systemLHS);
		}
		if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
			if (systemLHS instanceof SparseStore<?>) {
				int numberOfNonzero = Iterators.size(((SparseStore<?>) systemLHS).nonzeros());
				mainLog.println(String.format("Solving linear equation of size %dx%d, %d non-zero entries",
						systemLHS.countRows(), systemLHS.countColumns(), numberOfNonzero), PrismLog.VL_HIGH);
			} else {
				mainLog.println(String.format("Solving linear equation of size %dx%d", systemLHS.countRows(), systemLHS.countColumns()), PrismLog.VL_HIGH);
			}
		}
		try {
			MatrixStore<Double> solution = solver.solve(systemLHS, systemRHS);
			assert systemLHS.multiply(solution).subtract(systemRHS).aggregateAll(Aggregator.MAXIMUM) < PrismUtils.epsilonDouble;
			if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
				if (useSvd) {
					mainLog.println(String.format("Condition of the system was %g", ((SingularValue<?>) solver).getCondition()));
				} else {
					mainLog.println(String.format("Determinant of the system was %g", (Double) ((MatrixDecomposition.Determinant<?>) solver).getDeterminant()));
				}
			}
			return solution;
		} catch (TaskException e) {
			throw new PrismException("Could not solve linear equation system: " + e.getMessage());
		}

	}

	protected ProcessingResult<GainBiasResult> computeGainBiasLinearEquationsSccDecomposition(DTMC dtmc, MCRewards mcRewards, NatBitSet subset)
			throws PrismException
	{
		// SCC decomposition approach
		// Obtain gain and bias for each BSCC separately, then propagate backwards.

		Time time = new Time();

		if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
			mainLog.println("Computing gain and bias by linear equations and SCC decomposition", PrismLog.VL_HIGH);
		}

		int numDtmcStates = dtmc.getNumStates();
		int numStates = subset == null ? numDtmcStates : subset.size();

		ECComputerFast.DecompositionResult sccDecomposition = subset == null ? ECComputerFast.computeSCCs(dtmc)
				: ECComputerFast.computeSCCs(dtmc, subset.clone());

		assert sccDecomposition.sccs.stream().mapToInt(NatBitSet::size).sum() + sccDecomposition.bsccs.stream().mapToInt(NatBitSet::size).sum() == numStates;

		List<NatBitSet> sccs = sccDecomposition.sccs;
		List<NatBitSet> bsccs = sccDecomposition.bsccs;

		int numSccs = sccs.size();
		int numBsccs = bsccs.size();
		int numStatesInBsccs = bsccs.stream().mapToInt(NatBitSet::size).sum();
		int numStatesInSccs = numStates - numStatesInBsccs;
		boolean singleBscc = numSccs == 1;

		if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
			mainLog.println(String.format("Found %d SCCs (%d trivial, %d total states) and %d BSCCs (%d total states)", numSccs,
					sccs.stream().filter(scc -> scc.size() == 1).count(), numStatesInSccs, numBsccs, numStatesInBsccs), PrismLog.VL_HIGH);
		}

		Int2IntMap bsccStateToBsccIndexMap = new Int2IntLinkedOpenHashMap(numStatesInBsccs);
		bsccStateToBsccIndexMap.defaultReturnValue(-1);
		ListIterator<NatBitSet> bsccIterator = bsccs.listIterator();
		while (bsccIterator.hasNext()) {
			int bsccIndex = bsccIterator.nextIndex();
			NatBitSet bscc = bsccIterator.next();
			IntIterator bsccStateIterator = bscc.iterator();
			while (bsccStateIterator.hasNext()) {
				int bsccState = bsccStateIterator.nextInt();
				bsccStateToBsccIndexMap.put(bsccState, bsccIndex);
			}
		}

		Int2DoubleFunction biasMap;
		if (subset == null) {
			biasMap = new Nat2DoubleDenseArrayMap(numDtmcStates);
		} else {
			biasMap = new Int2DoubleLinkedOpenHashMap(numStates);
		}
		biasMap.defaultReturnValue(Double.NaN);
		double[] bsccGains = new double[numBsccs];
		Int2DoubleMap transientStateGainMap;
		if (singleBscc) {
			// All states have the same gain, no need to store it
			transientStateGainMap = null;
		} else {
			transientStateGainMap = new Int2DoubleLinkedOpenHashMap(numStatesInSccs);
			transientStateGainMap.defaultReturnValue(Double.NaN);
		}

		bsccIterator = bsccs.listIterator();
		while (bsccIterator.hasNext()) {
			int bsccIndex = bsccIterator.nextIndex();
			NatBitSet bscc = bsccIterator.next();
			GainBiasResult gainBias = computeGainBiasLinearEquationsUnichain(dtmc, mcRewards, bscc, bscc).getResult();

			bsccGains[bsccIndex] = gainBias.getGain(bscc.iterator().nextInt());
			assert bscc.stream().allMatch(s -> gainBias.getGain(s) == bsccGains[bsccIndex]);
			IntIterator bsccStateIterator = bscc.iterator();
			while (bsccStateIterator.hasNext()) {
				int bsccState = bsccStateIterator.nextInt();
				assert !biasMap.containsKey(bsccState) : bsccState;
				biasMap.put(bsccState, gainBias.getBias(bsccState));
			}
		}

		IntToDoubleFunction stateToGainFunction;
		if (singleBscc) {
			stateToGainFunction = state -> bsccGains[0];
		} else {
			stateToGainFunction = state -> {
				int bsccSuccessor = bsccStateToBsccIndexMap.get(state);
				if (bsccSuccessor != -1) {
					// This successor is in a bscc - take the gain of the bscc
					return bsccGains[bsccSuccessor];
				}
				// Successor is in a scc
				assert transientStateGainMap.containsKey(state);
				return transientStateGainMap.get(state);
			};
		}
		IntToDoubleFunction stateToBiasFunction = state -> {
			assert biasMap.containsKey(state);
			return biasMap.get(state);
		};

		NatBitSet bsccSuccessors = NatBitSets.boundedSet(numDtmcStates);
		NatBitSet sccSuccessorStates = NatBitSets.boundedSet(numDtmcStates);
		for (NatBitSet scc : sccs) {
			if (scc.size() == 1) {
				// Trivial SCC - No need for a equation system
				int sccState = scc.firstInt();

				// Gain: g(s) = sum_s' d(s, s') g(s') + d(s, s) g(s)
				// Bias: h(s) = sum_s' d(s, s') h(s') + d(s, s) h(s) + r(s) - g(s)

				double selfLoopProb = 0.0d;
				double gain = 0.0d;
				double bias = 0.0d;

				Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(sccState);
				while (transitionIterator.hasNext()) {
					Int2DoubleMap.Entry transition = transitionIterator.next();
					int destinationState = transition.getIntKey();

					double prob = transition.getDoubleValue();
					if (destinationState == sccState) {
						selfLoopProb += prob;
					} else {
						gain += prob * stateToGainFunction.applyAsDouble(destinationState);
						bias += prob * stateToBiasFunction.applyAsDouble(destinationState);
					}
				}

				bias += mcRewards.getStateReward(sccState) - gain;
				if (selfLoopProb > 0.0d) {
					assert selfLoopProb < 1.0d : "Not a SCC";
					gain /= 1 - selfLoopProb;
					bias /= 1 - selfLoopProb;
				}

				if (!singleBscc) {
					assert !transientStateGainMap.containsKey(sccState);
					transientStateGainMap.put(sccState, gain);
				}
				assert !biasMap.containsKey(sccState);
				biasMap.put(sccState, bias);
				continue;
			}

			IntUnaryOperator sccStateIndexMap = FastUtils.elementToIndexMap(scc);

			// Now iterate over all SCCs, successively back-propagating the values. As the scc list is topologically sorted, we have already determined
			// the values of all successor states.

			// First, solve reachability of possible successors for all states in the SCC by linear equations.
			// For this, we determine all successors and their gain values (states in other SCCs or BSCCs)

			IntIterator sccStateIterator = scc.iterator();
			while (sccStateIterator.hasNext()) {
				int sccState = sccStateIterator.nextInt();

				IntIterator successorsIterator = dtmc.getSuccessorsIterator(sccState);
				while (successorsIterator.hasNext()) {
					int successorState = successorsIterator.nextInt();

					if (scc.contains(successorState)) {
						continue;
					}
					int bsccSuccessor = bsccStateToBsccIndexMap.get(successorState);
					if (bsccSuccessor != -1) {
						bsccSuccessors.set(bsccSuccessor);
					} else {
						sccSuccessorStates.set(successorState);
					}
				}
			}
			IntUnaryOperator bsccSuccessorIndexMap = FastUtils.elementToIndexMap(bsccSuccessors);
			IntUnaryOperator sccSuccessorStateIndexMap = FastUtils.elementToIndexMap(sccSuccessorStates);
			IntList bsccSuccessorList = new IntArrayList(bsccSuccessors.size());
			IntList sccSuccessorStateList = new IntArrayList(sccSuccessorStates.size());
			bsccSuccessorList.addAll(bsccSuccessors);
			sccSuccessorStateList.addAll(sccSuccessorStates);
			// Cleanup
			bsccSuccessors.clear();
			sccSuccessorStates.clear();

			double[] bsccSuccessorGains = new double[bsccSuccessorList.size()];
			double[] sccSuccessorGains = new double[sccSuccessorStateList.size()];
			DoubleSet successorGains = new DoubleLinkedOpenHashSet();

			int numSccStates = scc.size();
			int numBsccSuccessors = bsccSuccessorList.size();
			int numSccSuccessorStates = sccSuccessorStateList.size();
			int numSuccessors = numBsccSuccessors + numSccSuccessorStates;
			assert numSuccessors > 0;

			IntListIterator bsccSuccessorIterator = bsccSuccessorList.listIterator();
			while (bsccSuccessorIterator.hasNext()) {
				int bsccSuccessorIndex = bsccSuccessorIterator.nextIndex();
				int bsccSuccessor = bsccSuccessorIterator.nextInt();
				double gain = bsccGains[bsccSuccessor];
				bsccSuccessorGains[bsccSuccessorIndex] = gain;
				successorGains.add(gain);
			}
			IntListIterator sccSuccessorIterator = sccSuccessorStateList.listIterator();
			while (sccSuccessorIterator.hasNext()) {
				int sccSuccessorStateIndex = sccSuccessorIterator.nextIndex();
				int sccSuccessorState = sccSuccessorIterator.nextInt();
				assert transientStateGainMap != null && transientStateGainMap.containsKey(sccSuccessorState);
				double gain = transientStateGainMap.get(sccSuccessorState);
				sccSuccessorGains[sccSuccessorStateIndex] = gain;
				successorGains.add(gain);
			}
			DoubleToIntFunction successorGainToIndexMap = Indices.elementToIndexMap(successorGains);

			int numDifferentGainSuccessors = successorGains.size();
			boolean trivialGain = numDifferentGainSuccessors == 1;
			assert !singleBscc || numSuccessors == 1;
			// If numSuccessors == 1 the gain must be trivial
			assert !(numSuccessors == 1) || trivialGain;

			// Now we create a RHS vector for each different successor gain.
			// Equation reach = P reach + reach_succ, i.e. (P - I) reach = - reach_gain
			SparseStore<Double> systemLHS = SparseStore2.makePrimitive(numSccStates, numSccStates);
			SparseStore<Double> systemRHS = trivialGain ? null : SparseStore2.makePrimitive(numSccStates, numDifferentGainSuccessors);

			sccStateIterator = scc.iterator();
			while (sccStateIterator.hasNext()) {
				int sccState = sccStateIterator.nextInt();
				int sccStateIndex = sccStateIndexMap.applyAsInt(sccState);

				Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(sccState);
				while (transitionIterator.hasNext()) {
					Int2DoubleMap.Entry transition = transitionIterator.next();
					int destinationState = transition.getIntKey();
					double prob = transition.getDoubleValue();

					if (scc.contains(destinationState)) {
						// P
						int destinationSccStateIndex = sccStateIndexMap.applyAsInt(destinationState);
						systemLHS.add(sccStateIndex, destinationSccStateIndex, prob);
					} else {
						if (trivialGain) {
							continue;
						}
						// -REACH
						int destinationIndex;
						// This is the "global" index of the bscc, i.e. its position in the bscc list
						int successorBsccIndex = bsccStateToBsccIndexMap.get(destinationState);
						if (successorBsccIndex != -1) {
							// This is the bscc's index in the rhs
							int bsccSuccessorIndex = bsccSuccessorIndexMap.applyAsInt(successorBsccIndex);
							double bsccGain = bsccSuccessorGains[bsccSuccessorIndex];
							destinationIndex = successorGainToIndexMap.applyAsInt(bsccGain);
						} else {
							assert sccSuccessorStateList.contains(destinationState);
							int sccSuccessorStateIndex = sccSuccessorStateIndexMap.applyAsInt(destinationState);
							double stateGain = sccSuccessorGains[sccSuccessorStateIndex];
							destinationIndex = successorGainToIndexMap.applyAsInt(stateGain);
						}
						assert destinationIndex != -1;
						systemRHS.add(sccStateIndex, destinationIndex, -prob);
					}
				}

				// -I
				systemLHS.add(sccStateIndex, sccStateIndex, -1.0d);
			}

			LU<Double> solver = LU.make(systemLHS);
			if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
				mainLog.println(String.format("Decomposing %dx%d transition matrix (%d non-zero entries)", numSccStates, numSccStates,
						Iterators.size(systemLHS.nonzeros())), PrismLog.VL_HIGH);
			}
			solver.decompose(systemLHS);

			if (!singleBscc) {
				if (trivialGain) {
					assert numDifferentGainSuccessors == 1;
					double gain = successorGains.iterator().nextDouble();

					sccStateIterator = scc.iterator();
					while (sccStateIterator.hasNext()) {
						int sccState = sccStateIterator.nextInt();
						assert !transientStateGainMap.containsKey(sccState);
						transientStateGainMap.put(sccState, gain);
					}
				} else {
					MatrixStore<Double> gainReachabilitySolution = solver.getSolution(systemRHS);

					// Fill the gain array with the obtained values
					sccStateIterator = scc.iterator();
					while (sccStateIterator.hasNext()) {
						int sccState = sccStateIterator.nextInt();
						int sccStateIndex = sccStateIndexMap.applyAsInt(sccState);

						double totalReachProb = 0.0d;
						double gain = 0.0d;
						DoubleIterator successorGainIterator = successorGains.iterator();
						while (successorGainIterator.hasNext()) {
							double successorGain = successorGainIterator.nextDouble();
							int successorGainIndex = successorGainToIndexMap.applyAsInt(successorGain);
							double gainReachability = gainReachabilitySolution.doubleValue(sccStateIndex, successorGainIndex);
							if (gainReachability == 0.0d) {
								continue;
							}
							gain += successorGain * gainReachability;
							totalReachProb += gainReachability;
						}
						assert PrismUtils.doublesAreEqual(totalReachProb, 1.0d);
						assert !transientStateGainMap.containsKey(sccState);
						transientStateGainMap.put(sccState, gain);
					}
				}
			}

			// Solve the bias equation system. We have that r - g + (P - I) h = 0 and that g(s) = sum_succ reach(succ) * gain(succ)
			// so we can move r - g to the RHS, overall its (P - I) h = - P * h_succ - r + g.

			// We can reuse the decomposition of the previous system, as it is P - I again so only build new RHS
			PrimitiveDenseStore systemBiasRHS = PrimitiveDenseStore.FACTORY.makeZero(numStates, 1);

			sccStateIterator = scc.iterator();
			while (sccStateIterator.hasNext()) {
				int sccState = sccStateIterator.nextInt();
				int sccStateIndex = sccStateIndexMap.applyAsInt(sccState);

				Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(sccState);
				while (transitionIterator.hasNext()) {
					Int2DoubleMap.Entry transition = transitionIterator.next();
					int destinationState = transition.getIntKey();
					double prob = transition.getDoubleValue();

					if (scc.contains(destinationState)) {
						continue;
					}
					double successorBias = stateToBiasFunction.applyAsDouble(destinationState);
					systemBiasRHS.add(sccStateIndex, -prob * successorBias);
				}

				double rhs = -mcRewards.getStateReward(sccState);
				if (singleBscc) {
					rhs += bsccGains[0];
				} else {
					assert transientStateGainMap.containsKey(sccState);
					rhs += transientStateGainMap.get(sccState);
				}
				systemBiasRHS.add(sccStateIndex, rhs);
			}

			MatrixStore<Double> biasValues = solver.getSolution(systemBiasRHS);

			sccStateIterator = scc.iterator();
			while (sccStateIterator.hasNext()) {
				int sccState = sccStateIterator.nextInt();
				int sccStateIndex = sccStateIndexMap.applyAsInt(sccState);

				assert !biasMap.containsKey(sccState);
				biasMap.put(sccState, biasValues.doubleValue(sccStateIndex));
			}
		}

		GainBiasResult gainBias = new GainBiasResultSparse(stateToGainFunction, stateToBiasFunction);
		return new ProcessingResult<>(gainBias, PrismUtils.epsilonDouble, 1, time.elapsedSeconds());
	}

	protected ProcessingResult<GainBiasResult> computeGainBiasLinearEquationsBsccCompression(DTMC dtmc, MCRewards mcRewards, NatBitSet subset)
			throws PrismException
	{
		Time time = new Time();

		if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
			mainLog.println("Computing gain and bias by BSCC compression", PrismLog.VL_HIGH);
		}

		int numDtmcStates = dtmc.getNumStates();
		int numStates = subset == null ? numDtmcStates : subset.size();

		List<NatBitSet> bsccs;
		if (subset == null) {
			bsccs = ECComputerFast.computeBSCCs(dtmc);
		} else {
			bsccs = ECComputerFast.computeBSCCs(dtmc, subset.clone());
		}

		if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
			mainLog.println(String.format("Found %d BSCCs, computing attractors", bsccs.size()), PrismLog.VL_HIGH);
		}

		NatBitSet transientStates;
		List<NatBitSet> attractors = new ArrayList<>(bsccs.size());

		if (bsccs.size() == 1) {
			transientStates = NatBitSets.emptySet();
			NatBitSet attractor = subset == null ? NatBitSets.fullSet(dtmc.getNumStates()) : subset.clone();
			attractors.add(attractor);
		} else {
			if (subset == null) {
				transientStates = NatBitSets.boundedFilledSet(numDtmcStates);
			} else {
				transientStates = subset.clone();
				assert Iterators.elementsEqual(transientStates.iterator(), subset.iterator());
			}
			attractors.addAll(bsccs);
			bsccs.forEach(transientStates::andNot);

			if (transientStates.isEmpty()) {
				if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
					mainLog.println("No transient states left (all states belong to a BSCC or some attractor)", PrismLog.VL_HIGH);
				}
			} else {
				PredecessorRelation predecessorRelation;
				if (subset == null) {
					predecessorRelation = dtmc.getPredecessorRelation(this, false);
				} else {
					predecessorRelation = new PredecessorRelationSparse(dtmc, subset);
				}
				ListIterator<NatBitSet> iterator = bsccs.listIterator();
				while (iterator.hasNext()) {
					int bsccIndex = iterator.nextIndex();
					NatBitSet bscc = iterator.next();
					NatBitSet bsccAttractor = prob1(dtmc, transientStates, bscc, predecessorRelation);
					if (subset != null) {
						bsccAttractor.and(subset);
					}
					assert bsccAttractor.containsAll(bscc);
					if (bscc.size() == bsccAttractor.size()) {
						continue;
					}
					attractors.set(bsccIndex, bsccAttractor);
					transientStates.andNot(bsccAttractor);

					if (transientStates.isEmpty()) {
						if (mainLog.isLoggable(PrismLog.VL_HIGH)) {
							mainLog.println("No transient states left (all states belong to a BSCC or some attractor)", PrismLog.VL_HIGH);
						}
						break;
					}
				}
			}
		}

		IntUnaryOperator transientStateIndexMap = FastUtils.elementToIndexMap(transientStates);
		int numTransientStates = transientStates.size();
		int numStatesInAttractor = numStates - numTransientStates;

		int stateIndex = 0;
		Int2IntMap attractorStateToBsccIndexMap = new Int2IntLinkedOpenHashMap(numStatesInAttractor);
		attractorStateToBsccIndexMap.defaultReturnValue(Integer.MIN_VALUE);
		Int2IntMap attractorStateToStateIndexMap = new Int2IntLinkedOpenHashMap(numStatesInAttractor);
		attractorStateToStateIndexMap.defaultReturnValue(Integer.MIN_VALUE);
		ListIterator<NatBitSet> attractorIterator = attractors.listIterator();
		while (attractorIterator.hasNext()) {
			int attractorIndex = attractorIterator.nextIndex();
			NatBitSet attractor = attractorIterator.next();
			IntIterator attractorStateIterator = attractor.iterator();
			while (attractorStateIterator.hasNext()) {
				int attractorState = attractorStateIterator.nextInt();
				assert !attractorStateToBsccIndexMap.containsKey(attractorState);
				attractorStateToBsccIndexMap.put(attractorState, attractorIndex);
				assert !attractorStateToStateIndexMap.containsKey(attractorState);
				attractorStateToStateIndexMap.put(attractorState, stateIndex);
				stateIndex += 1;
			}
		}
		assert stateIndex == numStatesInAttractor;

		assert subset == null || subset.containsAll(transientStates);
		assert subset == null || attractors.stream().allMatch(subset::containsAll);

		/* Approach here: We compress each BSCC (and its attractor) into a single gain variable and sort them on top. So
		 * the gain vector will be (bscc_1, bscc_2, ..., bscc_n, trans_1, ..., trans_m), where trans are transient states.
		 * The bias vector instead will contain one entry for each state, but with corresponding ordering, i.e.
		 * (bscc_1_1, ..., bscc_1_k1, ..., bscc_2_1, ..., bscc_n_kn, trans_1, ..., trans_m).
		 *
		 * Since the top left entries for the collapsed BSCC states cancel out (they only have self loops)
		 * we don't need to consider them.
		 */

		int numBsccs = bsccs.size();
		// One gain variable per transient state and bscc, one bias variable per state
		int gainRows = numTransientStates;
		int gainColumns = numBsccs + numTransientStates;

		// One row for each gain variable, one row for each bias variable, and one row to fix the bias variables
		int rows = gainRows + numStates + numBsccs;
		int columns = gainColumns + numStates;
		// Columns are: <gain attr><gain transient><bias attr><bias transient>
		// Width:        numBssc    numTransient    numAttr    numTransient
		// Rows are: <gain transient><bias attr><bias transient>
		// Height:    numTransient    numAttr    numTransient
		SparseStore<Double> systemLHS = SparseStore2.makePrimitive(rows, columns);
		PrimitiveDenseStore systemRHS = PrimitiveDenseStore.FACTORY.makeZero(rows, 1);

		// Set the transient stuff
		IntIterator transientStateIterator = transientStates.iterator();
		while (transientStateIterator.hasNext()) {
			int transientState = transientStateIterator.nextInt();
			int transientStateIndex = transientStateIndexMap.applyAsInt(transientState);

			int transientStateGainRow = transientStateIndex;
			int transientStateBiasRow = gainRows + numStatesInAttractor + transientStateIndex;
			int transientStateGainColumn = numBsccs + transientStateIndex;
			int transientStateBiasColumn = gainColumns + numStatesInAttractor + transientStateIndex;

			Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(transientState);
			while (transitionIterator.hasNext()) {
				Int2DoubleMap.Entry transition = transitionIterator.next();
				int destinationState = transition.getIntKey();
				assert subset == null || subset.contains(destinationState);
				double prob = transition.getDoubleValue();
				if (transientStates.contains(destinationState)) {
					int destinationIndex = transientStateIndexMap.applyAsInt(destinationState);
					// Top left P
					systemLHS.set(transientStateGainRow, numBsccs + destinationIndex, prob);
					// Bottom right P
					systemLHS.set(transientStateBiasRow, gainColumns + numStatesInAttractor + destinationIndex, prob);
				} else {
					// Top left P
					assert attractorStateToBsccIndexMap.containsKey(destinationState);
					int destinationAttractorIndex = attractorStateToBsccIndexMap.get(destinationState);
					systemLHS.add(transientStateGainRow, destinationAttractorIndex, prob);

					// Bottom right P
					assert attractorStateToStateIndexMap.containsKey(destinationState);
					int destinationStateIndex = attractorStateToStateIndexMap.get(destinationState);
					systemLHS.set(transientStateBiasRow, gainColumns + destinationStateIndex, prob);
				}
			}
			// Add top left -I
			systemLHS.add(transientStateGainRow, transientStateGainColumn, -1.0d);
			// Add bottom right -I
			systemLHS.add(transientStateBiasRow, transientStateBiasColumn, -1.0d);
			// Set bottom left -gain
			systemLHS.set(transientStateBiasRow, transientStateGainColumn, -1.0d);

			// Set rhs
			systemRHS.set(transientStateBiasRow, -mcRewards.getStateReward(transientState));
		}

		// Set bottom left -gain and bottom right bias-P for all attractor states
		attractorIterator = attractors.listIterator();
		while (attractorIterator.hasNext()) {
			int attractorIndex = attractorIterator.nextIndex();
			NatBitSet attractor = attractorIterator.next();
			IntIterator attractorStateIterator = attractor.iterator();
			while (attractorStateIterator.hasNext()) {
				int attractorState = attractorStateIterator.nextInt();
				int attractorStateIndex = attractorStateToStateIndexMap.get(attractorState);

				int attractorStateGainColumn = attractorIndex;
				int attractorStateBiasRow = gainRows + attractorStateIndex;
				int attractorStateBiasColumn = gainColumns + attractorStateIndex;

				Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(attractorState);
				while (transitionIterator.hasNext()) {
					Int2DoubleMap.Entry transition = transitionIterator.next();
					int destinationState = transition.getIntKey();
					double prob = transition.getDoubleValue();

					assert attractor.contains(destinationState) && attractorStateToStateIndexMap.containsKey(destinationState);
					int destinationIndex = attractorStateToStateIndexMap.get(destinationState);
					systemLHS.set(attractorStateBiasRow, gainColumns + destinationIndex, prob);
				}
				// Bottom left -I
				systemLHS.set(attractorStateBiasRow, attractorStateGainColumn, -1.0d);
				// Bottom right -I
				systemLHS.add(attractorStateBiasRow, attractorStateBiasColumn, -1.0d);
				// RHS
				systemRHS.set(attractorStateBiasRow, -mcRewards.getStateReward(attractorState));
			}
		}

		// Get unique solution by adding bias constraints
		ListIterator<NatBitSet> bsccIterator = bsccs.listIterator();
		while (bsccIterator.hasNext()) {
			int bsccIndex = bsccIterator.nextIndex();
			NatBitSet bscc = bsccIterator.next();
			int bsccBiasFixedState = bscc.iterator().nextInt();
			int bsccBiasFixedStateIndex = attractorStateToStateIndexMap.get(bsccBiasFixedState);
			int row = gainRows + numStates + bsccIndex;
			systemLHS.set(row, gainColumns + bsccBiasFixedStateIndex, 1.0d);
			systemRHS.set(row, 0.0d);
		}

		MatrixStore<Double> solution = solveLinearEquationSystem(systemLHS, systemRHS, false);
		double[] gainArray = new double[gainColumns];
		double[] biasArray = new double[numStates];

		assert gainArray.length + biasArray.length == Iterables.size(solution.rows());
		for (int i = 0; i < gainColumns; i++) {
			gainArray[i] = solution.doubleValue(i);
		}
		for (int i = 0; i < numStates; i++) {
			biasArray[i] = solution.doubleValue(gainColumns + i);
		}

		IntToDoubleFunction stateToGainFunction = state -> {
			if (transientStates.contains(state)) {
				return gainArray[numBsccs + transientStateIndexMap.applyAsInt(state)];
			}
			return gainArray[attractorStateToBsccIndexMap.get(state)];
		};
		IntToDoubleFunction stateToBiasFunction = state -> {
			if (transientStates.contains(state)) {
				return biasArray[numStatesInAttractor + transientStateIndexMap.applyAsInt(state)];
			}
			return biasArray[attractorStateToStateIndexMap.get(state)];
		};

		GainBiasResult gainBias = new GainBiasResultSparse(stateToGainFunction, stateToBiasFunction);
		return new ProcessingResult<>(gainBias, PrismUtils.epsilonDouble, 1, time.elapsedSeconds());
	}

	protected ProcessingResult<GainBiasResult> computeGainBiasLinearEquations(DTMC dtmc, MCRewards mcRewards, NatBitSet subset)
			throws PrismException
	{
		// TODO Make the lin eq library plugable

		if (EnumSet.of(DTMC.Type.UNICHAIN, DTMC.Type.STRONGLY_UNICHAIN).contains(dtmc.getMCType())) {
			return computeGainBiasLinearEquationsUnichain(dtmc, mcRewards, subset, null);
		}

		Time time = new Time();

		List<NatBitSet> bsccs = subset == null ? ECComputerFast.computeBSCCs(dtmc) : ECComputerFast.computeBSCCs(dtmc, subset.clone());

		int numDtmcStates = dtmc.getNumStates();
		int numStates = subset == null ? numDtmcStates : subset.size();
		int numBsccs = bsccs.size();

		// Solve (P - I)g = 0 along with r - g + (P - I)h = 0 and h(s) = 0 for a state s for each BSCC
		// | P-I(1) |    0(4) | | g |    0(7)
		// |  -I(2) |  P-I(5) | | h | = -r(8)
		// |   0(3) | BSCC(6) |          0(9)
		int rowCount = numStates * 2 + numBsccs;

		SparseStore<Double> systemLHS = SparseStore2.makePrimitive(rowCount, numStates * 2);
		PrimitiveDenseStore systemRHS = PrimitiveDenseStore.FACTORY.makeZero(rowCount, 1);

		// Build the equation system row by row.
		IntUnaryOperator stateToIndexMap = FastUtils.elementToIndexMap(subset);
		IntIterator stateIterator = FastUtils.iterator(subset, numDtmcStates);
		while (stateIterator.hasNext()) {
			int state = stateIterator.nextInt();
			int stateIndex = stateToIndexMap.applyAsInt(state);

			Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(state);
			while (transitionIterator.hasNext()) {
				Int2DoubleMap.Entry transition = transitionIterator.next();
				int destinationIndex = stateToIndexMap.applyAsInt(transition.getIntKey());
				double prob = transition.getDoubleValue();

				if (prob == 0.0d) {
					continue;
				}
				// (1) P
				systemLHS.set(stateIndex, destinationIndex, prob);
				// (5) P
				systemLHS.set(numStates + stateIndex, numStates + destinationIndex, prob);
			}
			// (1) -I
			systemLHS.add(stateIndex, stateIndex, -1.0d);
			// (2)
			systemLHS.set(numStates + stateIndex, stateIndex, -1.0d);
			// (5) -I
			systemLHS.add(numStates + stateIndex, numStates + stateIndex, -1.0d);

			// (8)
			double reward = mcRewards.getStateReward(state);
			if (reward != 0) {
				systemRHS.set(numStates + stateIndex, -reward);
			}
		}

		int zeroStateNumber = 0;
		for (NatBitSet bscc : bsccs) {
			int zeroState = bscc.firstInt();
			int zeroStateIndex = stateToIndexMap.applyAsInt(zeroState);
			// (6)
			systemLHS.set(2L * numStates + zeroStateNumber, numStates + zeroStateIndex, 1.0d);
			zeroStateNumber += 1;
		}

		MatrixStore<Double> solution = solveLinearEquationSystem(systemLHS, systemRHS, true);
		double[] gainArray = new double[numStates];
		double[] biasArray = new double[numStates];
		for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
			// First numStates rows correspond to gain and the remaining correspond to bias
			double gainValue = solution.doubleValue(stateIndex);
			double biasValue = solution.doubleValue(numStates + stateIndex);
			assert !Double.isNaN(gainValue) && !Double.isNaN(biasValue);
			gainArray[stateIndex] = gainValue;
			biasArray[stateIndex] = biasValue;
		}

		GainBiasResult gainBias;
		if (subset == null) {
			gainBias = new GainBiasResultDense(gainArray, biasArray);
		} else {
			gainBias = new GainBiasResultSparse(
					state -> gainArray[stateToIndexMap.applyAsInt(state)],
					state -> biasArray[stateToIndexMap.applyAsInt(state)]);
		}
		return new ProcessingResult<>(gainBias, PrismUtils.epsilonDouble, 1, time.elapsedSeconds());
	}

	protected ProcessingResult<GainBiasResult> computeGainBiasLinearEquationsUnichain(DTMC dtmc, MCRewards mcRewards, NatBitSet subset, NatBitSet bscc)
			throws PrismException
	{
		Time time = new Time();

		int numDtmcStates = dtmc.getNumStates();

		if (bscc == null) {
			List<NatBitSet> bsccs = ECComputerFast.computeBSCCs(dtmc, subset);
			checkArgument(bsccs.size() == 1);
			bscc = bsccs.get(0);
		}

		int numStates = subset == null ? numDtmcStates : subset.size();
		int numBsccStates = bscc.size();
		int rowCount = numBsccStates + 1;

		// Solve r - g1 + (P - I)h = 0 on the BSCC
		SparseStore<Double> systemLHS = SparseStore2.makePrimitive(rowCount, numBsccStates + 1);
		PrimitiveDenseStore systemRHS = PrimitiveDenseStore.FACTORY.makeZero(rowCount, 1);
		IntUnaryOperator bsccStateToIndexMap = FastUtils.elementToIndexMap(bscc);

		IntIterator bsccIterator = bscc.iterator();
		while (bsccIterator.hasNext()) {
			int bsccState = bsccIterator.nextInt();
			int bsccStateIndex = bsccStateToIndexMap.applyAsInt(bsccState);

			Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(bsccState);
			// -g
			systemLHS.set(bsccStateIndex, 0, -1.0d);
			while (transitionIterator.hasNext()) {
				Int2DoubleMap.Entry transition = transitionIterator.next();
				int destinationState = transition.getIntKey();
				double prob = transition.getDoubleValue();
				assert bscc.contains(destinationState);

				int destinationStateIndex = bsccStateToIndexMap.applyAsInt(destinationState);
				// P for bias
				systemLHS.set(bsccStateIndex, 1 + destinationStateIndex, prob);
			}
			// -I for bias
			systemLHS.add(bsccStateIndex, 1 + bsccStateIndex, -1.0d);
			// -r RHS
			systemRHS.set(bsccStateIndex, -mcRewards.getStateReward(bsccState));
		}

		int zeroState = bscc.firstInt();
		int zeroStateIndex = bsccStateToIndexMap.applyAsInt(zeroState);
		// Uniqueness for bias
		systemLHS.set(numBsccStates, 1 + zeroStateIndex, 1);

		MatrixStore<Double> bsccSolution = solveLinearEquationSystem(systemLHS, systemRHS, false);
		double bsccGain = bsccSolution.doubleValue(0);
		double[] bsccBiasArray = new double[numBsccStates];
		Arrays.setAll(bsccBiasArray, s -> bsccSolution.doubleValue(1 + s));

		if (numBsccStates == numStates) {
			// We are done here
			GainBiasResult gainBias;
			if (subset == null) {
				gainBias = new GainBiasResultSingle(bsccGain, bsccBiasArray);
			} else {
				gainBias = new GainBiasResultSparse(state -> bsccGain, state -> bsccBiasArray[bsccStateToIndexMap.applyAsInt(state)]);
			}
			return new ProcessingResult<>(gainBias, PrismUtils.epsilonDouble, 1, time.elapsedSeconds());
		}

		// Now build the system for all transient states and solve for bias (we know the gain already)
		// Solve (P_trans - I) h = -P_bscc * h_bscc + g1 - r
		NatBitSet transientStates;
		if (subset == null) {
			transientStates = NatBitSets.boundedFilledSet(numStates);
		} else {
			transientStates = subset.clone();
		}
		transientStates.andNot(bscc);

		int numTransientStates = transientStates.size();
		IntUnaryOperator transientStateToIndexMap = FastUtils.elementToIndexMap(transientStates);
		systemLHS = SparseStore2.makePrimitive(numTransientStates, numTransientStates);
		systemRHS = PrimitiveDenseStore.FACTORY.makeZero(numTransientStates, 1);

		IntIterator transientStateIterator = transientStates.iterator();

		while (transientStateIterator.hasNext()) {
			int transientState = transientStateIterator.nextInt();
			int transientStateIndex = transientStateToIndexMap.applyAsInt(transientState);

			Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(transientState);
			while (transitionIterator.hasNext()) {
				Int2DoubleMap.Entry transition = transitionIterator.next();
				int destinationState = transition.getIntKey();
				double prob = transition.getDoubleValue();

				if (bscc.contains(destinationState)) {
					// RHS -P_bscc * h_bscc
					systemRHS.add(transientStateIndex, -prob * bsccBiasArray[bsccStateToIndexMap.applyAsInt(destinationState)]);
				} else {
					int destinationStateIndex = transientStateToIndexMap.applyAsInt(destinationState);
					systemLHS.set(transientStateIndex, destinationStateIndex, prob);
				}
			}
			// -I for bias
			systemLHS.add(transientStateIndex, transientStateIndex, -1.0d);
			// Rhs gain - r
			systemRHS.add(transientStateIndex, bsccGain - mcRewards.getStateReward(transientState));
		}

		MatrixStore<Double> solution = solveLinearEquationSystem(systemLHS, systemRHS, false);

		double[] transientBiasArray = new double[numTransientStates];
		Arrays.setAll(transientBiasArray, solution::doubleValue);

		NatBitSet finalBscc = bscc;
		GainBiasResult gainBias = new GainBiasResultSparse(state -> bsccGain, state -> {
			if (finalBscc.contains(state)) {
				return bsccBiasArray[bsccStateToIndexMap.applyAsInt(state)];
			}
			return transientBiasArray[transientStateToIndexMap.applyAsInt(state)];
		});
		return new ProcessingResult<>(gainBias, PrismUtils.epsilonDouble, 1, time.elapsedSeconds());
	}

	/**
	 * Compute transient probabilities
	 * i.e. compute the probability of being in each state at time step {@code k},
	 * assuming the initial distribution {@code initDist}.
	 * For space efficiency, the initial distribution vector will be modified and values over-written,
	 * so if you wanted it, take a copy.
	 *
	 * @param dtmc     The DTMC
	 * @param k        Time step
	 * @param initDist Initial distribution (will be overwritten)
	 */
	public ModelCheckerResult computeTransientProbs(DTMC dtmc, int k, double initDist[]) throws PrismException
	{
		throw new PrismNotSupportedException("Not implemented yet");
	}

	public ModelCheckerResult computeAverageReward(DTMC dtmc, MCRewards modelRewards) throws PrismException
	{
		double precision = termCritParam;

		if (settings.getString(PrismSettings.PRISM_MC_MP_METHOD).equals("Interval iteration")) {
			return IntervalResult.toModelCheckerResult(computeAverageRewardIntervalValueIteration(dtmc, modelRewards, null, precision), dtmc.getNumStates());
		} else if (settings.getString(PrismSettings.PRISM_MC_MP_METHOD).equals("Linear equations")) {
			return GainResult.toModelCheckerResult(computeGainBias(dtmc, modelRewards, null), dtmc.getNumStates());
		} else {
			throw new PrismException("Solution method " + settings.getString(PrismSettings.PRISM_MC_MP_METHOD) + " not supported");
		}
	}

	ProcessingResult<IntervalResult> computeAverageRewardIntervalValueIteration(DTMC dtmc, MCRewards modelRewards, NatBitSet subset, double precision)
			throws PrismException
	{
		Time time = new Time();
		MCSpanNormValueIterator iterator = new MCSpanNormValueIterator(this, (DTMCExplicit) dtmc, modelRewards);
		if (dtmc.getMCType() == STRONGLY_UNICHAIN || dtmc.getMCType() == UNICHAIN) {
			MCSpanNormValueIterator.MCRunner runner = iterator.createRunner(null, precision);
			runner.run(-1);
			return new ProcessingResult<>(new IntervalResultSingle(runner.getResultMinimumValue(), runner.getResultMaximumValue()),
					runner.getObtainedPrecision(), runner.getIterations(), runner.getTimeTaken());
		}

		// Find BSCCs
		List<NatBitSet> bsccs = ECComputerFast.computeBSCCs(dtmc, subset);
		double minimalGain = Double.POSITIVE_INFINITY;
		double maximalGain = Double.NEGATIVE_INFINITY;
		double[] bsccMinimalGains = new double[bsccs.size()];
		double[] bsccMaximalGains = new double[bsccs.size()];
		int valueIterations = 0;

		// Evaluate for each bscc
		int bsccIndex = 0;
		for (NatBitSet bscc : bsccs) {
			MCSpanNormValueIterator.MCRunner runner = iterator.createRunner(bscc, precision / 2);
			runner.run(-1);
			valueIterations += runner.getIterations();
			double bsccMaximalGain = runner.getResultMaximumValue();
			double bsccMinimalGain = runner.getResultMinimumValue();
			assert bsccMinimalGain <= bsccMaximalGain;
			bsccMinimalGains[bsccIndex] = bsccMinimalGain;
			bsccMaximalGains[bsccIndex] = bsccMaximalGain;
			if (bsccMinimalGain < minimalGain) {
				minimalGain = bsccMinimalGain;
			}
			if (bsccMaximalGain > maximalGain) {
				maximalGain = bsccMaximalGain;
			}
			bsccIndex += 1;
		}

		double gainDifference = maximalGain - minimalGain;
		assert gainDifference >= 0;
		if (gainDifference < precision) {
			// Every BSCC yields practically the same bound
			IntervalResultSingle intervalResult = new IntervalResultSingle(minimalGain, maximalGain);
			return new ProcessingResult<>(intervalResult, gainDifference, valueIterations, time.elapsedSeconds());
		}

		// Collapse the DTMC
		DTMCCollapsed collapsedDtmc = new DTMCCollapsed(dtmc, subset, bsccs);

		int numCollapsedStates = collapsedDtmc.getNumStates();
		double[] reachLower = new double[numCollapsedStates];
		double[] reachUpper = new double[numCollapsedStates];
		double[] reachLowerPrevious = new double[numCollapsedStates];
		double[] reachUpperPrevious = new double[numCollapsedStates];
		Arrays.fill(reachUpper, 1d);
		Arrays.fill(reachUpperPrevious, 1d);

		NatBitSet unknown = NatBitSets.boundedFilledSet(numCollapsedStates);
		bsccIndex = 0;
		for (NatBitSet bscc : bsccs) {
			int firstBsccState = bscc.iterator().nextInt();
			int bsccCollapsedState = collapsedDtmc.getCollapsedStateIndex(firstBsccState);
			reachLower[bsccCollapsedState] = (bsccMinimalGains[bsccIndex] - minimalGain) / gainDifference;
			reachUpper[bsccCollapsedState] = (bsccMaximalGains[bsccIndex] - minimalGain) / gainDifference;
			unknown.clear(bsccCollapsedState);
			bsccIndex += 1;
		}
		// TODO prob0 step here? compute the prob1 of all bsccs with minimalGain (+- epsilon) and mark them as cleared
		double reachabilityPrecision = precision / 2 * gainDifference;

		// Start iterations
		int reachabilityIterations = 0;
		boolean done = false;
		while (!done && reachabilityIterations < maxIters) {
			reachabilityIterations++;

			collapsedDtmc.mvMult(reachLower, reachLowerPrevious, unknown);
			double[] swapLower = reachLower;
			reachLower = reachLowerPrevious;
			reachLowerPrevious = swapLower;
			collapsedDtmc.mvMult(reachUpper, reachUpperPrevious, unknown);
			double[] swapUpper = reachUpper;
			reachUpper = reachUpperPrevious;
			reachUpperPrevious = swapUpper;

			done = true;
			for (int state = 0; state < reachUpper.length; state++) {
				if (!PrismUtils.doublesAreCloseAbs(reachLower[state], reachUpper[state], reachabilityPrecision)) {
					done = false;
					break;
				} else {
					unknown.clear(state);
				}
			}
		}

		// Convert it back into a vector of values
		double finalMinimalGain = minimalGain;
		double[] finalReachUpper = reachUpper;
		IntToDoubleFunction upperFunction = state -> {
			int collapsedState = collapsedDtmc.getCollapsedStateIndex(state);
			return finalReachUpper[collapsedState] * gainDifference + finalMinimalGain;
		};
		double[] finalReachLower = reachLower;
		IntToDoubleFunction lowerFunction = state -> {
			int collapsedState = collapsedDtmc.getCollapsedStateIndex(state);
			return finalReachLower[collapsedState] * gainDifference;
		};
		assert Iterators.all(IntIterators.fromTo(0, numCollapsedStates), state -> finalReachLower[state] <= finalReachUpper[state]);

		IntervalResultSparse intervalResult = new IntervalResultSparse(lowerFunction, upperFunction);
		return new ProcessingResult<>(intervalResult, precision, valueIterations + reachabilityIterations, time.elapsedSeconds());
	}

	protected IntervalValueIterator createGainValueIterator(DTMC model, MCRewards rewards, IntIterable subset)
	{
		boolean sparse = subset != null;
		double maxReward = 0.0d;
		IntIterator subsetIterator = subset == null ? IntIterators.fromTo(0, model.getNumStates()) : subset.iterator();
		while (subsetIterator.hasNext()) {
			int nextState = subsetIterator.nextInt();
			double reward = rewards.getStateReward(nextState);
			if (reward > maxReward) {
				maxReward = reward;
			}
		}

		if (sparse) {
			return new GainValueIterator(model, rewards, maxReward, Int2DoubleOpenHashMap::new);
		}
		int numStates = model.getNumStates();
		Supplier<Int2DoubleFunction> storageCreator = () -> new Nat2DoubleDenseArrayMap(numStates);
		return new GainValueIterator(model, rewards, maxReward, storageCreator);
	}

	public ProcessingResult<GainBiasResult> computeGainBias(DTMC dtmc, MCRewards mcRewards, NatBitSet subset) throws PrismException
	{
		if (dtmc.getMCType() == STRONGLY_UNICHAIN) {
			return computeGainBiasLinearEquationsUnichain(dtmc, mcRewards, subset, subset);
		}
		if (dtmc.getMCType() == UNICHAIN) {
			return computeGainBiasLinearEquationsUnichain(dtmc, mcRewards, subset, null);
		}

		String solvingMethodString = settings.getString(PrismSettings.PRISM_MC_GAIN_BIAS_METHOD);
		AverageRewardLinearEquationMethod solvingMethod;
		if (solvingMethodString.equals("Full")) {
			solvingMethod = AverageRewardLinearEquationMethod.FULL_SVD;
		} else if (solvingMethodString.equals("BSCC compression")) {
			solvingMethod = AverageRewardLinearEquationMethod.BSCC_COMPRESSION;
		} else if (solvingMethodString.equals("SCC decomposition")) {
			solvingMethod = AverageRewardLinearEquationMethod.SCC_DECOMPOSITION;
		} else {
			throw new PrismException("Unknown solving method " + solvingMethodString);
		}

		switch (solvingMethod) {
		case FULL_SVD:
			return computeGainBiasLinearEquations(dtmc, mcRewards, subset);
		case BSCC_COMPRESSION:
			return computeGainBiasLinearEquationsBsccCompression(dtmc, mcRewards, subset);
		case SCC_DECOMPOSITION:
			return computeGainBiasLinearEquationsSccDecomposition(dtmc, mcRewards, subset);
		default:
			throw new PrismException("Unknown solving method " + solvingMethod);
		}
	}

	public enum AverageRewardEvaluationMethod
	{
		LINEAR_EQUATION, INTERVAL_ITERATION, VALUE_ITERATION
	}

	public enum AverageRewardSolutionMethod
	{
		FULL, BSCC_DECOMPOSITION
	}

	public enum AverageRewardDecompositionReachabilityMethod
	{
		VALUE_ITERATION, INTERVAL_ITERATION
	}

	public enum AverageRewardLinearEquationMethod
	{
		FULL_SVD, BSCC_COMPRESSION, SCC_DECOMPOSITION
	}

	interface IntervalValueIterator
	{
		void setLowerBound(int state, double value);

		void setUpperBound(int state, double value);

		double getLowerBound(int state);

		double getUpperBound(int state);

		boolean run(NatBitSet subset, int steps, double precision);
	}

	static class GainValueIterator implements IntervalValueIterator
	{
		private final DTMC model;
		private final MCRewards rewards;
		private Int2DoubleFunction currentLowerValues;
		private Int2DoubleFunction currentUpperValues;
		private Int2DoubleFunction previousLowerValues;
		private Int2DoubleFunction previousUpperValues;

		GainValueIterator(DTMC model, MCRewards rewards, double upperBound, Supplier<Int2DoubleFunction> storageSupplier)
		{
			this.model = model;
			this.rewards = rewards;
			this.currentLowerValues = storageSupplier.get();
			this.currentUpperValues = storageSupplier.get();
			this.previousLowerValues = storageSupplier.get();
			this.previousUpperValues = storageSupplier.get();
			currentLowerValues.defaultReturnValue(0);
			previousLowerValues.defaultReturnValue(0);
			currentUpperValues.defaultReturnValue(upperBound);
			previousUpperValues.defaultReturnValue(upperBound);
		}

		@Override public void setLowerBound(int state, double value)
		{
			currentLowerValues.put(state, value);
		}

		@Override public void setUpperBound(int state, double value)
		{
			currentUpperValues.put(state, value);
		}

		@Override public double getLowerBound(int state)
		{
			return currentLowerValues.get(state) - previousLowerValues.get(state);
		}

		@Override public double getUpperBound(int state)
		{
			return currentUpperValues.get(state) - previousUpperValues.get(state);
		}

		private void swap()
		{
			Int2DoubleFunction swap = currentLowerValues;
			currentLowerValues = previousLowerValues;
			previousLowerValues = swap;
			swap = currentUpperValues;
			currentUpperValues = previousUpperValues;
			previousUpperValues = swap;
		}

		@Override public boolean run(NatBitSet subset, int stepBound, double precision)
		{
			int numDtmcStates = model.getNumStates();
			int numStates = subset == null ? numDtmcStates : subset.size();

			NatBitSet statesToProcess = subset == null ? NatBitSets.boundedFilledSet(numDtmcStates) : subset.clone();

			IntIterator iterator = FastUtils.iterator(subset, numDtmcStates);

			int steps = 0;
			while (!statesToProcess.isEmpty() && steps < stepBound) {
				IntIterator stateIterator = statesToProcess.iterator();

				swap();
				while (stateIterator.hasNext()) {
					int state = stateIterator.nextInt();

					double stateReward = rewards.getStateReward(state);
					double previousLowerValue = previousLowerValues.get(state);
					double previousUpperValue = previousUpperValues.get(state);
					double currentLowerValue = stateReward + model.mvMultSingle(state, previousLowerValues::get);
					double currentUpperValue = stateReward + model.mvMultSingle(state, previousUpperValues::get);
					currentLowerValues.put(state, currentLowerValue);
					currentUpperValues.put(state, currentUpperValue);

					if (steps == 0) {
						continue;
					}

					double upperDifference = currentUpperValue - previousUpperValue;
					double lowerDifference = currentLowerValue - previousLowerValue;

					assert PrismUtils.doublesAreLessOrEqual(0, lowerDifference);
					assert PrismUtils.doublesAreLessOrEqual(0, upperDifference);
					assert PrismUtils.doublesAreLessOrEqual(lowerDifference, upperDifference);
					if (currentUpperValue - currentLowerValue <= precision) {
						// statesToProcess.clear(state);
					}
				}
				steps++;
			}
			return statesToProcess.isEmpty();
		}
	}

	static class MCSpanNormValueIterator extends explicit.SpanNormValueIterator
	{
		private final DTMCExplicit model;
		private final MCRewards rewards;

		MCSpanNormValueIterator(PrismComponent parent, DTMCExplicit model, MCRewards rewards)
		{
			super(parent);
			this.model = model;
			this.rewards = rewards;
		}

		@Override Model getModel()
		{
			return model;
		}

		public MCRunner createRunner(NatBitSet scc, double precision)
		{
			return new MCRunner(scc, step -> step % 5 == 1, precision);
		}

		public MCRunner createRunner(NatBitSet scc, IntPredicate checkStepPredicate, double precision)
		{
			return new MCRunner(scc, checkStepPredicate, precision);
		}

		class MCRunner extends Runner
		{
			private final NatBitSet scc;

			MCRunner(NatBitSet scc, IntPredicate checkStepPredicate, double precision)
			{
				super(checkStepPredicate, precision);
				this.scc = scc;
			}

			@Override double operator(int state, IntToDoubleFunction stateValues, double tau)
			{
				assert model.allSuccessorsInSet(state, scc);
				return tau * rewards.getStateReward(state) + model.mvMultSingleAperiodic(state, stateValues, tau);
			}

			@Override IntIterator getStateIterator()
			{
				return FastUtils.iterator(scc, model.getNumStates());
			}
		}
	}
}
