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
import com.google.common.collect.Iterators;
import common.FastUtils;
import common.Time;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.result.GainBiasResult;
import explicit.result.GainResult;
import explicit.result.GainResultSingle;
import explicit.result.GainResultSparse;
import explicit.result.IntervalResult;
import explicit.result.ProcessingResult;
import explicit.rewards.MCRewards;
import explicit.rewards.MCRewardsFromMDPRewards;
import explicit.rewards.MDPRewards;
import explicit.rewards.Rewards;
import it.unimi.dsi.fastutil.doubles.DoubleComparator;
import it.unimi.dsi.fastutil.doubles.DoubleComparators;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectFunction;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntLists;
import it.unimi.dsi.fastutil.ints.IntSet;
import parser.VarList;
import parser.ast.Declaration;
import parser.ast.DeclarationIntUnbounded;
import parser.ast.Expression;
import prism.Prism;
import prism.PrismComponent;
import prism.PrismDevNullLog;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import prism.PrismSettings;
import prism.PrismUtils;
import strat.MDStrategy;
import strat.MDStrategyArray;
import strat.MDStrategySparse;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.function.IntConsumer;
import java.util.function.IntPredicate;
import java.util.function.IntToDoubleFunction;
import java.util.function.Supplier;

/**
 * Explicit-state model checker for Markov decision processes (MDPs).
 */
public class MDPModelChecker extends ProbModelChecker
{
	/**
	 * Create a new MDPModelChecker, inherit basic state from parent (unless null).
	 */
	public MDPModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	// Model checking functions

	/**
	 * Simple test program.
	 */
	public static void main(String args[])
	{
		MDPModelChecker mc;
		MDPSimple mdp;
		ModelCheckerResult res;
		NatBitSet init, target;
		Map<String, NatBitSet> labels;
		boolean min = true;
		try {
			mc = new MDPModelChecker(null);
			mdp = new MDPSimple();
			mdp.buildFromPrismExplicit(args[0]);
			mdp.addInitialState(0);
			//System.out.println(mdp);
			labels = StateModelChecker.loadLabelsFile(args[1]);
			//System.out.println(labels);
			init = labels.get("init");
			target = labels.get(args[2]);
			if (target == null)
				throw new PrismException("Unknown label \"" + args[2] + "\"");
			for (int i = 3; i < args.length; i++) {
				if (args[i].equals("-min"))
					min = true;
				else if (args[i].equals("-max"))
					min = false;
				else if (args[i].equals("-nopre"))
					mc.setPrecomp(false);
			}
			res = mc.computeReachProbs(mdp, target, min);
			System.out.println(res.soln[init.firstInt()]);
		} catch (PrismException e) {
			System.out.println(e);
		}
	}

	@Override
	protected StateValues checkProbPathFormulaLTL(Model model, Expression expr, boolean qual, MinMax minMax, NatBitSet statesOfInterest) throws PrismException
	{
		StateValues probsProduct, probs;

		// For min probabilities, need to negate the formula
		// (add parentheses to allow re-parsing if required)
		if (minMax.isMin()) {
			expr = Expression.Not(Expression.Parenth(expr.deepCopy()));
		}

		boolean ldba = settings.getString(PrismSettings.PRISM_MDP_LTL_METHOD).equals("LDBA");

		Product<MDP> product;
		// For LTL model checking routines
		LTLModelChecker modelChecker = new LTLModelChecker(this);

		if (ldba) {
			product = modelChecker.constructProductMDPLDBA(this, (MDP) model, expr, statesOfInterest);
		} else {
			// Build product of MDP and automaton
			AcceptanceType[] allowedAcceptance = {
					AcceptanceType.BUCHI,
					AcceptanceType.RABIN,
					AcceptanceType.GENERALIZED_RABIN,
					AcceptanceType.REACH
			};
			product = modelChecker.constructProductMDP(this, (MDP) model, expr, statesOfInterest, allowedAcceptance);
		}

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

		if (ldba) {
			// Find accepting states + compute reachability probabilities
			mainLog.println("\nComputing reachability probabilities...");
			MDPModelChecker mcProduct = new MDPModelChecker(this);
			mcProduct.inheritSettings(this);
			probsProduct = StateValues.createFromDoubleArray(mcProduct.computeReachProbs(product.getProductModel(),
					((LTLModelChecker.LDBAProduct<?>) product).getGoalStates(), false).soln, product.getProductModel());
		} else {
			LTLModelChecker.LTLProduct<?> ltlProduct = (LTLModelChecker.LTLProduct<?>) product;

			// Find accepting states + compute reachability probabilities
			NatBitSet acc;
			if (ltlProduct.getAcceptance() instanceof AcceptanceReach) {
				mainLog.println("\nSkipping accepting MEC computation since acceptance is defined via goal states...");
				acc = ((AcceptanceReach) ltlProduct.getAcceptance()).getGoalStates();
			} else {
				mainLog.println("\nFinding accepting MECs...");
				acc = modelChecker.findAcceptingECStates(product.getProductModel(), ltlProduct.getAcceptance());
			}
			mainLog.println("\nComputing reachability probabilities...");
			MDPModelChecker mcProduct = new MDPModelChecker(this);
			mcProduct.inheritSettings(this);
			ModelCheckerResult res = mcProduct.computeReachProbs(product.getProductModel(), acc, false);
			probsProduct = StateValues.createFromDoubleArray(res.soln, product.getProductModel());
		}

		// Subtract from 1 if we're model checking a negated formula for regular Pmin
		if (minMax.isMin()) {
			probsProduct.timesConstant(-1.0);
			probsProduct.plusConstant(1.0);
		}

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

	// Numerical computation functions

	/**
	 * Compute rewards for a co-safe LTL reward operator.
	 */
	@Override protected StateValues checkRewardCoSafeLTL(Model model, Rewards modelRewards, Expression expr, MinMax minMax, NatBitSet statesOfInterest)
			throws PrismException
	{
		LTLModelChecker mcLtl;
		MDPRewards productRewards;
		StateValues rewardsProduct, rewards;
		MDPModelChecker mcProduct;
		LTLModelChecker.LTLProduct<MDP> product;

		// For LTL model checking routines
		mcLtl = new LTLModelChecker(this);

		// Build product of MDP and automaton
		AcceptanceType[] allowedAcceptance = {
				AcceptanceType.RABIN,
				AcceptanceType.REACH
		};
		product = mcLtl.constructProductMDP(this, (MDP) model, expr, statesOfInterest, allowedAcceptance);

		// Adapt reward info to product model
		productRewards = ((MDPRewards) modelRewards).liftFromModel(product);

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
			// For a DFA, just collect the accept states
			mainLog.println("\nSkipping end component detection since DRA is a DFA...");
			acc = ((AcceptanceReach) product.getAcceptance()).getGoalStates();
		} else {
			// Usually, we have to detect end components in the product
			mainLog.println("\nFinding accepting end components...");
			acc = mcLtl.findAcceptingECStates(product.getProductModel(), product.getAcceptance());
		}
		mainLog.println("\nComputing reachability rewards...");
		mcProduct = new MDPModelChecker(this);
		mcProduct.inheritSettings(this);
		ModelCheckerResult res = mcProduct.computeReachRewards(product.getProductModel(), productRewards, acc, minMax.isMin());
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

	/**
	 * Compute next=state probabilities.
	 * i.e. compute the probability of being in a state in {@code target} in the next step.
	 *
	 * @param mdp    The MDP
	 * @param target Target states
	 * @param min    Min or max probabilities (true=min, false=max)
	 */
	public ModelCheckerResult computeNextProbs(MDP mdp, NatBitSet target, boolean min) throws PrismException
	{
		ModelCheckerResult res;
		int n;
		double soln[], soln2[];
		long timer;

		timer = System.currentTimeMillis();

		// Store num states
		n = mdp.getNumStates();

		// Create/initialise solution vector(s)
		soln = Utils.intSetToDoubleArray(target, n);
		soln2 = new double[n];

		// Next-step probabilities
		mdp.mvMultMinMax(soln, min, soln2, NatBitSets.fullSet(mdp.getNumStates()), null);

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln2;
		res.numIters = 1;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Given a value vector x, compute the probability:
	 * v(s) = min/max sched [ Sum_s' P_sched(s,s')*x(s') ]  for s labeled with a,
	 * v(s) = 0   for s not labeled with a.
	 * <p>
	 * Clears the StateValues object x.
	 *
	 * @param mdp the transition matrix
	 * @param a   the set of states labeled with a
	 * @param x   the value vector
	 * @param min compute min instead of max
	 */
	public double[] computeRestrictedNext(MDP mdp, NatBitSet a, double[] x, boolean min)
	{
		int n;
		double soln[];

		// Store num states
		n = mdp.getNumStates();

		// initialized to 0.0
		soln = new double[n];

		// Next-step probabilities multiplication
		// restricted to a states
		mdp.mvMultMinMax(x, min, soln, a, null);

		return soln;
	}

	/**
	 * Compute reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target}.
	 *
	 * @param mdp    The MDP
	 * @param target Target states
	 * @param min    Min or max probabilities (true=min, false=max)
	 */
	public ModelCheckerResult computeReachProbs(MDP mdp, NatBitSet target, boolean min) throws PrismException
	{
		return computeReachProbs(mdp, null, target, min, null, null);
	}

	/**
	 * Compute until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in {@code remain}.
	 *
	 * @param mdp    The MDP
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min    Min or max probabilities (true=min, false=max)
	 */
	public ModelCheckerResult computeUntilProbs(MDP mdp, NatBitSet remain, NatBitSet target, boolean min) throws PrismException
	{
		return computeReachProbs(mdp, remain, target, min, null, null);
	}

	public ModelCheckerResult computeReachProbs(MDP mdp, NatBitSet remain, NatBitSet target, boolean min, double init[], IntSet known, MDPSolnMethod method)
			throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet no, yes;
		int n, numYes, numNo;
		long timer, timerProb0, timerProb1;
		MDStrategy strat = null;
		// Local copy of setting
		MDPSolnMethod mdpSolnMethod = method;

		// Check for some unsupported combinations
		if (mdpSolnMethod == MDPSolnMethod.VALUE_ITERATION && valIterDir == ValIterDir.ABOVE) {
			if (!(precomp && prob0))
				throw new PrismException("Precomputation (Prob0) must be enabled for value iteration from above");
			if (!min)
				throw new PrismException("Value iteration from above only works for minimum probabilities");
		}
		if (mdpSolnMethod == MDPSolnMethod.POLICY_ITERATION || mdpSolnMethod == MDPSolnMethod.MODIFIED_POLICY_ITERATION) {
			if (known != null) {
				throw new PrismException("Policy iteration methods cannot be passed 'known' values for some states");
			}
		}

		// Start probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting probabilistic reachability (" + (min ? "min" : "max") + ")...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		mdp.checkForDeadlocks(target);

		// Store num states
		n = mdp.getNumStates();

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null && !known.isEmpty()) {
			NatBitSet targetNew = target.clone();
			known.forEach((IntConsumer) i -> {
				if (init[i] == 1.0) {
					targetNew.add(i);
				}
			});
			target = targetNew;
		}

		// If required, export info about target states
		if (getExportTarget()) {
			NatBitSet bsInit = NatBitSets.boundedSet(n);
			for (int i = 0; i < n; i++) {
				bsInit.set(i, mdp.isInitialState(i));
			}
			List<NatBitSet> labels = Arrays.asList(bsInit, target);
			List<String> labelNames = Arrays.asList("init", "target");
			mainLog.println("\nExporting target states info to file \"" + getExportTargetFilename() + "\"...");
			exportLabels(mdp, labels, labelNames, Prism.EXPORT_PLAIN, new PrismFileLog(getExportTargetFilename()));
		}

		// If required, create/initialise strategy storage
		// Set choices to -1, denoting unknown
		// (except for target states, which are -2, denoting arbitrary)
		if (genStrat || exportAdv) {
			strat = new MDStrategyArray(mdp, new int[n]);
			for (int i = 0; i < n; i++) {
				strat.setChoiceIndex(i, target.contains(i) ? -2 : -1);
			}
		}

		// Precomputation
		timerProb0 = System.currentTimeMillis();
		if (precomp && prob0) {
			no = prob0(mdp, remain, target, min, strat);
		} else {
			no = NatBitSets.emptySet();
		}
		timerProb0 = System.currentTimeMillis() - timerProb0;
		timerProb1 = System.currentTimeMillis();
		if (precomp && prob1) {
			yes = prob1(mdp, remain, target, min, strat);
		} else {
			yes = target.clone();
		}
		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numYes = yes.size();
		numNo = no.size();
		mainLog.println("target=" + target.size() + ", yes=" + numYes + ", no=" + numNo + ", maybe=" + (n - (numYes + numNo)));

		// If still required, store strategy for no/yes (0/1) states.
		// This is just for the cases max=0 and min=1, where arbitrary choices suffice (denoted by -2)
		if (genStrat || exportAdv) {
			assert strat != null;
			if (min) {
				IntIterator iterator = yes.iterator();
				while (iterator.hasNext()) {
					int i = iterator.nextInt();
					if (!target.contains(i)) {
						strat.setChoiceIndex(i, -2);
					}
				}
			} else {
				IntIterator iterator = no.iterator();
				while (iterator.hasNext()) {
					int i = iterator.nextInt();
					strat.setChoiceIndex(i, -2);
				}
			}
		}

		// Compute probabilities (if needed)
		if (numYes + numNo < n) {
			switch (mdpSolnMethod) {
			case VALUE_ITERATION:
				res = computeReachProbsValIter(mdp, no, yes, min, init, known, strat);
				break;
			case GAUSS_SEIDEL:
				res = computeReachProbsGaussSeidel(mdp, no, yes, min, init, known, strat);
				break;
			case POLICY_ITERATION:
				res = computeReachProbsPolIter(mdp, no, yes, min, strat);
				break;
			case MODIFIED_POLICY_ITERATION:
				res = computeReachProbsModPolIter(mdp, no, yes, min, strat);
				break;
			case INTERVAL_ITERATION:
				res = computeReachProbsIntervalIter(mdp, no, yes, min, strat);
				break;
			default:
				throw new PrismException("Unknown MDP solution method " + mdpSolnMethod.fullName());
			}
		} else {
			res = new ModelCheckerResult();
			res.soln = Utils.intSetToDoubleArray(yes, n);
		}

		// Finished probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Probabilistic reachability took " + timer / 1000.0 + " seconds.");

		// Store strategy
		if (genStrat) {
			res.strat = strat;
		}
		// Export adversary
		if (exportAdv) {
			// Prune strategy
			restrictStrategyToReachableStates(mdp, strat);
			// Export
			PrismLog out = new PrismFileLog(exportAdvFilename);
			new DTMCFromMDPAndMDStrategy(mdp, strat).exportToPrismExplicitTra(out);
			out.close();
		}

		// Update time taken
		res.timeTaken = timer / 1000.0;
		res.timeProb0 = timerProb0 / 1000.0;
		res.timePre = (timerProb0 + timerProb1) / 1000.0;

		return res;
	}

	/**
	 * Compute reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in {@code remain}.
	 *
	 * @param mdp    The MDP
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min    Min or max probabilities (true=min, false=max)
	 * @param init   Optionally, an initial solution vector (may be overwritten)
	 * @param known  Optionally, a set of states for which the exact answer is known
	 *               Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values).
	 *               Also, 'known' values cannot be passed for some solution methods, e.g. policy iteration.
	 */
	public ModelCheckerResult computeReachProbs(MDP mdp, NatBitSet remain, NatBitSet target, boolean min, double init[], NatBitSet known) throws PrismException
	{
		return computeReachProbs(mdp, remain, target, min, init, known, mdpSolnMethod);
	}

	/**
	 * Prob0 precomputation algorithm.
	 * i.e. determine the states of an MDP which, with min/max probability 0,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 * {@code min}=true gives Prob0E, {@code min}=false gives Prob0A.
	 * Optionally, for min only, store optimal (memoryless) strategy info for 0 states.
	 *
	 * @param mdp    The MDP
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min    Min or max probabilities (true=min, false=max)
	 * @param strat  Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	public NatBitSet prob0(MDP mdp, NatBitSet remain, NatBitSet target, boolean min, MDStrategy strat)
	{
		int n, iters;
		NatBitSet u, soln, unknown;
		boolean u_done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		mainLog.println("Starting Prob0 (" + (min ? "min" : "max") + ")...");

		// Special case: no target states
		if (target.isEmpty()) {
			return NatBitSets.fullSet(mdp.getNumStates());
		}

		// Initialise vectors
		n = mdp.getNumStates();
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
			mdp.prob0step(unknown, u, min, soln);
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
		mainLog.print("Prob0 (" + (min ? "min" : "max") + ")");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		// If required, generate strategy. This is for min probs,
		// so it can be done *after* the main prob0 algorithm (unlike for prob1).
		// We simply pick, for all "no" states, the first choice for which all transitions stay in "no"
		if (strat != null) {
			u.forEach((IntConsumer) i -> {
				int numChoices = mdp.getNumChoices(i);
				for (int k = 0; k < numChoices; k++) {
					if (mdp.allSuccessorsInSet(i, k, u)) {
						strat.setChoiceIndex(i, k);
						break;
					}
				}
			});
		}

		return u;
	}

	/**
	 * Prob1 precomputation algorithm.
	 * i.e. determine the states of an MDP which, with min/max probability 1,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 * {@code min}=true gives Prob1A, {@code min}=false gives Prob1E.
	 * Optionally, for max only, store optimal (memoryless) strategy info for 1 states.
	 *
	 * @param mdp    The MDP
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min    Min or max probabilities (true=min, false=max)
	 * @param strat  Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	public NatBitSet prob1(MDP mdp, NatBitSet remain, NatBitSet target, boolean min, MDStrategy strat)
	{
		// Start precomputation
		long timer = System.currentTimeMillis();
		mainLog.println("Starting Prob1 (" + (min ? "min" : "max") + ")...");

		// Special case: no target states
		if (target.isEmpty()) {
			return NatBitSets.boundedSet(mdp.getNumStates());
		}

		// Initialise vectors
		int n = mdp.getNumStates();
		// Greatest fixed point
		NatBitSet u = NatBitSets.boundedFilledSet(n);
		NatBitSet v = NatBitSets.boundedSet(n);
		NatBitSet soln = NatBitSets.boundedSet(n);

		// Determine set of states actually need to perform computation for
		NatBitSet unknown = NatBitSets.boundedFilledSet(n);
		unknown.andNot(target);
		if (remain != null) {
			unknown.and(remain);
		}

		// Nested fixed point loop
		int iters = 0;
		boolean u_done = false;
		boolean v_done;
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
				if (min) {
					mdp.prob1Astep(unknown, u, v, soln);
				} else {
					mdp.prob1Estep(unknown, u, v, soln, null);
				}
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

		// If we need to generate a strategy, do another iteration of the inner loop for this
		// We could do this during the main double fixed point above, but we would generate surplus
		// strategy info for non-1 states during early iterations of the outer loop,
		// which are not straightforward to remove since this method does not know which states
		// already have valid strategy info from Prob0.
		// Notice that we only need to look at states in u (since we already know the answer),
		// so we restrict 'unknown' further
		unknown.and(u);
		if (!min && strat != null) {
			v_done = false;
			v.clear();
			v.or(target);
			soln.clear();
			soln.or(target);
			while (!v_done) {
				mdp.prob1Estep(unknown, u, v, soln, strat);
				v_done = soln.equals(v);
				v.clear();
				v.or(soln);
			}
			u_done = v.equals(u);
		}

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Prob1 (" + (min ? "min" : "max") + ")");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		return u;
	}

	/**
	 * Compute reachability probabilities using value iteration.
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param mdp   The MDP
	 * @param no    Probability 0 states
	 * @param yes   Probability 1 states
	 * @param min   Min or max probabilities (true=min, false=max)
	 * @param init  Optionally, an initial solution vector (will be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 *              Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachProbsValIter(MDP mdp, IntSet no, IntSet yes, boolean min, double init[], IntSet known, MDStrategy strat)
			throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[], initVal;
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting value iteration (" + (min ? "min" : "max") + ")...");

		// Store num states
		n = mdp.getNumStates();

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
			// Matrix-vector multiply and min/max ops
			mdp.mvMultMinMax(soln, min, soln2, unknown, strat);
			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Value iteration (" + (min ? "min" : "max") + ")");
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
	 * Interval iteration as in
	 * Serge Haddad, Benjamin Monmege, Interval iteration algorithm for MDPs and IMDPs, Theoretical Computer Science,
	 * Available online 2 January 2017, ISSN 0304-3975, http://dx.doi.org/10.1016/j.tcs.2016.12.003.
	 * <p>
	 * Maintains upper (initialized to 1) and lower bounds (intialized to 0) for all states and performs
	 * value iteration on them. The update is performed both over the upper as well as the lower values.
	 * The difference between the upper and lower values gives the error bound.
	 *
	 * @param mdp             An MDP
	 * @param knownZeroStates Probability 0 states
	 * @param knownOneStates  Probability 1 states
	 * @param min             true for computing minimum reachability probability
	 * @param strat           Storage for memoryless deterministic strategy
	 * @return
	 */
	protected ModelCheckerResult computeReachProbsIntervalIter(MDP mdp, IntSet knownZeroStates, IntSet knownOneStates, boolean min, MDStrategy strat)
			throws PrismException
	{
		// Start value iteration
		Time time = new Time();
		mainLog.println(String.format("%nStarting reachability interval iteration (%s)...", min ? "min" : "max"));

		// Store num states
		int numOriginalStates = mdp.getNumStates();
		int numCollapsedStates;
		MDP collapsedModel;

		double[] previousLowerValues;
		double[] previousUpperValues;
		double[] currentLowerValues;
		double[] currentUpperValues;

		if (mdp.getMDPInformation().contains(MDP.Info.NO_TRANSIENT_MEC)) {
			// No need to collapse

			// Create solution vector(s)
			previousLowerValues = new double[numOriginalStates];
			previousUpperValues = new double[numOriginalStates];
			currentLowerValues = new double[numOriginalStates];
			currentUpperValues = new double[numOriginalStates];

			// Initialize
			Arrays.fill(previousUpperValues, 1d);
			IntIterator oneIterator = knownOneStates.iterator();
			while (oneIterator.hasNext()) {
				previousLowerValues[oneIterator.nextInt()] = 1d;
			}
			IntIterator zeroIterator = knownZeroStates.iterator();
			while (zeroIterator.hasNext()) {
				previousUpperValues[zeroIterator.nextInt()] = 0d;
			}
			collapsedModel = mdp;
			numCollapsedStates = numOriginalStates;
		} else {
			// Find MECs in the MDP
			mainLog.println("Discovering MECs... ");
			List<MEC> mecs = ECComputerFast.computeMECs(mdp);

			// Filter out trivial (with no action) MECs and
			// set mecReward = 1 if MEC contains a goal state;
			// else set mecReward = 0
			Int2DoubleMap mecRewards = new Int2DoubleArrayMap();
			Collection<MEC> nonTrivialMecs = new ArrayList<>();

			for (MEC mec : mecs) {
				if (mec.actions.isEmpty()) {
					continue;
				}
				nonTrivialMecs.add(mec);
				int state = mec.states.firstInt();
				if (mec.states.intersects(knownOneStates)) {
					mecRewards.put(state, 1d);
				} else {
					mecRewards.put(state, 0d);
				}
				// TODO: Can we use knownZeroStates to our advantage anywhere?
			}

			mainLog.println(String.format("Found %d MECs, %d non-trivial and %d goal MECs", mecs.size(), nonTrivialMecs.size(), mecRewards.size()));

			mainLog.println("Collapsing MECs in MDP...");
			// Construct the MEC-collapsed MDP and transform it with added '+' and '-' states so that reachability can be run on it.
			MDPCollapsed collapsedMdp = new MDPCollapsed(mdp, mecs, mecRewards, knownOneStates);
			collapsedMdp.findDeadlocks(true);

			numCollapsedStates = collapsedMdp.getNumStates();
			collapsedModel = collapsedMdp;

			// Create solution vector(s)
			previousLowerValues = new double[numCollapsedStates];
			previousUpperValues = new double[numCollapsedStates];
			currentLowerValues = new double[numCollapsedStates];
			currentUpperValues = new double[numCollapsedStates];

			// Initialize
			Arrays.fill(previousUpperValues, 1d);
			previousLowerValues[collapsedMdp.getGoalState()] = 1d;
			previousUpperValues[collapsedMdp.getTrapState()] = 0d;
		}

		double precision = termCritParam;
		double timePre = time.elapsedSeconds(true);

		// Start iterations
		int iterations = 0;
		boolean done = false;
		while (!done) {
			iterations++;
			// Matrix-vector multiply and min/max ops
			// strat stores the min/max action chosen when updating the lower value
			collapsedModel.mvMultMinMax(previousLowerValues, min, currentLowerValues, NatBitSets.fullSet(numCollapsedStates), strat);
			collapsedModel.mvMultMinMax(previousUpperValues, min, currentUpperValues, NatBitSets.fullSet(numCollapsedStates), strat);

			// Check termination
			done = PrismUtils.doublesAreClose(currentUpperValues, currentLowerValues, precision, false);

			// Swap vectors for next iter
			double[] swapLower = previousLowerValues;
			previousLowerValues = currentLowerValues;
			currentLowerValues = swapLower;

			double[] swapUpper = previousUpperValues;
			previousUpperValues = currentUpperValues;
			currentUpperValues = swapUpper;
		}

		// Set solution = average of lower and upper values
		double[] reachabilityProbability = new double[numOriginalStates];
		for (int originalState = 0; originalState < numOriginalStates; originalState++) {
			int collapsedState = collapsedModel == mdp ? originalState :
					((MDPCollapsed) collapsedModel).getCollapsedStateIndexForOriginalState(originalState);
			reachabilityProbability[originalState] = (previousLowerValues[collapsedState] + previousUpperValues[collapsedState]) / 2;
		}

		// Finished value iteration
		double timeTaken = time.elapsedSeconds();
		mainLog.print("Interval iteration (" + (min ? "min" : "max") + ")");
		mainLog.println(String.format(" took %d iterations and %f seconds.", iterations, timeTaken));

		// Return results
		ModelCheckerResult result = new ModelCheckerResult();
		result.soln = reachabilityProbability;
		result.numIters = iterations;
		result.timeTaken = timeTaken;
		result.timePre = timePre;
		return result;
	}

	public ModelCheckerResult computeAverageReward(MDP mdp, MDPRewards rewards, boolean min) throws PrismException
	{
		ModelCheckerResult res;
		String solutionMethodString = settings.getString(PrismSettings.PRISM_MDP_SOLN_METHOD);
		if (solutionMethodString.equals("Value iteration")) {
			mainLog.println("Starting value iteration");
			mainLog.printWarning("Full VI has no convergence guarantees for general MDPs and should be used with great care");
			res = computeAverageRewardVI(mdp, rewards, min, null, null);
		} else if (solutionMethodString.equals("Policy iteration")) {
			mainLog.println("Starting policy iteration");
			ProcessingResult<? extends GainResult> piResult = computeAverageRewardPI(mdp, rewards, min, null, null);
			return GainResult.toModelCheckerResult(piResult, mdp.getNumStates());
		} else if (solutionMethodString.equals("MEC decomposition")) {
			AverageRewardMECSolutionMethod solutionMethod;
			String mecSolutionMethodString = settings.getString(PrismSettings.PRISM_MDP_MP_MEC_SOLN_METHOD);
			if (mecSolutionMethodString.equals("Value iteration")) {
				solutionMethod = AverageRewardMECSolutionMethod.VALUE_ITERATION;
			} else if (mecSolutionMethodString.equals("Policy iteration")) {
				solutionMethod = AverageRewardMECSolutionMethod.POLICY_ITERATION;
			} else {
				throw new PrismException("Unknown MEC solution method " + mecSolutionMethodString);
			}

			AverageRewardMECReachMethod reachMethod;
			String reachMethodString = settings.getString(PrismSettings.PRISM_MDP_MEC_REACH_METHOD);
			if (reachMethodString.equals("Interval iteration")) {
				reachMethod = AverageRewardMECReachMethod.INTERVAL_VALUE_ITERATION;
			} else if (reachMethodString.equals("Linear programming")) {
				reachMethod = AverageRewardMECReachMethod.LINEAR_PROGRAM;
			} else if (reachMethodString.equals("Value iteration")) {
				reachMethod = AverageRewardMECReachMethod.VALUE_ITERATION;
			} else if (reachMethodString.equals("Policy iteration")) {
				reachMethod = AverageRewardMECReachMethod.POLICY_ITERATION;
			} else {
				throw new PrismException("Unknown MEC reachability method " + reachMethodString);
			}

			mainLog.println(String.format("Starting MEC decomposition MP solver with mec solver %s and reach solver %s", solutionMethod, reachMethod));
			res = computeAverageRewardMECDecomposition(mdp, rewards, solutionMethod, reachMethod, min);
		} else {
			throw new PrismException("Method " + solutionMethodString + " not supported for mean payoff...");
		}
		return res;
	}

	/**
	 * Compute expected average rewards for a mec
	 *
	 * @param mdpRewards    The rewards
	 * @param min           Min or max rewards (true=min, false=max)
	 * @param mec           The MEC
	 * @param initialValues An initial vector of values for mec (or any subset of the mec)
	 */
	public ValueIterationModelCheckerResult computeAverageRewardVI(MDP mdp, MDPRewards mdpRewards, boolean min, MEC mec, double[] initialValues)
			throws PrismException
	{
		double precision = termCritParam;
		MDPSpanNormValueIterator valueIterator =
				new MDPSpanNormValueIterator(this, (MDPExplicit) mdp, mdpRewards, min);
		MDPSpanNormValueIterator.MDPRunner runner = valueIterator.createRunner(mec, precision);
		runner.run(-1);

		ValueIterationModelCheckerResult result = new ValueIterationModelCheckerResult();
		result.soln = new double[mdp.getNumStates()];
		FastUtils.iterator(mec, mdp.getNumStates()).forEachRemaining((IntConsumer) s -> result.soln[s] = runner.getCurrentGain(s));
		result.numIters = runner.getIterations();
		result.timeTaken = runner.getTimeTaken();
		result.meanPayoff = runner.getResultMeanPayoff();
		result.maximumValue = runner.getResultMaximumValue();
		result.minimumValue = runner.getResultMinimumValue();
		mainLog.println(String.format("Finished value iteration on MEC in %d iterations, result %f and %f secs.",
				result.numIters, result.meanPayoff, result.timeTaken));
		return result;
	}

	/**
	 * Compute expected average (step-bounded) rewards by decomposing the MDP into MECs, computing the
	 * gain for each MEC and then using weighted reachability to find the expected maximum average reward.
	 *
	 * @param mdpRewards The rewards
	 * @param min        Min or max rewards (true=min, false=max)
	 */
	public ModelCheckerResult computeAverageRewardMECDecomposition(MDP model, MDPRewards mdpRewards, AverageRewardMECSolutionMethod mecSolutionMethod,
			AverageRewardMECReachMethod reachabilitySolutionMethod, boolean min) throws PrismException
	{
		if (model.getNumInitialStates() > 1) {
			throw new UnsupportedOperationException("Can't run MP with multiple initial states.");
		}

		ModelCheckerResult res = new ModelCheckerResult();
		Time time = new Time();
		Time preTime = new Time();

		// Find MECs in the MDP
		mainLog.println("Discovering MECs... ");
		List<MEC> mecs = ECComputerFast.computeMECs(model);
		mainLog.println(String.format("Found %d mecs after %s seconds", mecs.size(), PrismUtils.formatDouble(time.elapsedSeconds())), PrismLog.VL_HIGH);
		double maxMECReward = 0d;
		Int2DoubleMap mecRewards = new Int2DoubleArrayMap();
		Collection<MEC> nonTrivialMecs = new ArrayList<>();

		// For each MEC, run the solver until convergence is reached.
		for (MEC mec : mecs) {
			if (mec.size() > 1) {
				nonTrivialMecs.add(mec);
				continue;
			}
			if (mec.actions.isEmpty()) {
				continue;
			}
			int state = mec.states.firstInt();
			assert mec.actions.size() == 1 && mec.actions.containsKey(state);
			nonTrivialMecs.add(mec);

			double maximalSelfLoopReward = Double.NEGATIVE_INFINITY;
			IntIterator actions = mec.actions.get(state).iterator();
			while (actions.hasNext()) {
				int action = actions.nextInt();
				maximalSelfLoopReward = Math.max(maximalSelfLoopReward, mdpRewards.getTransitionReward(state, action));
			}
			assert maximalSelfLoopReward >= 0d;
			maximalSelfLoopReward += mdpRewards.getStateReward(state);
			mecRewards.put(state, maximalSelfLoopReward);
			maxMECReward = Math.max(maxMECReward, maximalSelfLoopReward);
			res.numIters += 1;
		}
		res.timePre = preTime.elapsedSeconds(true);

		mainLog.println(String.format("Found %d MECs, %d non-trivial and pre-computed %d payoffs", mecs.size(), nonTrivialMecs.size(), mecRewards.size()));

		double precision = termCritParam;

		MecAverageRewardSolverFactory factory;
		if (mecSolutionMethod == AverageRewardMECSolutionMethod.VALUE_ITERATION) {
			factory = new ValueIterationSolverFactory();
		} else if (mecSolutionMethod == AverageRewardMECSolutionMethod.POLICY_ITERATION) {
			factory = new PolicyIterationSolverFactory();
		} else {
			throw new PrismException("Unsupported solution method " + mecSolutionMethod);
		}

		MecAverageRewardSolver solver = factory.obtain(model, mdpRewards, settings, min, precision);
		for (MEC mec : nonTrivialMecs) {
			int firstMECState = mec.states.firstInt();
			if (mecRewards.containsKey(firstMECState)) {
				// Processed previously when identified as a non-trivial single-state MEC
				continue;
			}
			mainLog.println(String.format("Processing MEC of size %d", mec.size()), PrismLog.VL_HIGH);
			MeanPayoffModelCheckerResult result = solver.solve(mec);
			res.numIters += result.numIters;
			mecRewards.put(firstMECState, result.meanPayoff);
			maxMECReward = Math.max(maxMECReward, result.meanPayoff);
		}

		res.soln = new double[model.getNumStates()];

		if (nonTrivialMecs.size() == 1) {
			// There only is a single MEC and every run will eventually end up in it. Hence we can immediately return the value
			mainLog.println(String.format("Returning value %f for single MEC", maxMECReward));
			Arrays.fill(res.soln, maxMECReward);
			res.timeTaken = time.elapsedSeconds();
			return res;
		}
		if (maxMECReward < precision / 2d) {
			// By assumption all rewards are non-negative, hence all MEC rewards are in [0, eps / 2) and the MP will be small.
			mainLog.println("Returning zero value");
			Arrays.fill(res.soln, maxMECReward / 2d);
			res.timeTaken = time.elapsedSeconds();
			return res;
		}

		preTime.restart();
		mainLog.println(String.format("Determined max reward %f, collapsing MECs in MDP and transforming the mean-payoff objective into a "
				+ "reachability objective...", maxMECReward));
		// Construct the MEC-collapsed MDP and transform it with added '+' and '-' states so that reachability can be run on it.
		MDPCollapsed collapsedMdp = new MDPCollapsed(model, mecs, mecRewards);
		collapsedMdp.findDeadlocks(true);

		NatBitSet goal = NatBitSets.boundedSet(collapsedMdp.getNumStates());
		// Set the '+' state as goal state
		goal.set(collapsedMdp.getGoalState());

		MDPSolnMethod reachabilityMethod;
		switch (reachabilitySolutionMethod) {
		case INTERVAL_VALUE_ITERATION:
			reachabilityMethod = MDPSolnMethod.INTERVAL_ITERATION;
			break;
		case LINEAR_PROGRAM:
			reachabilityMethod = MDPSolnMethod.LINEAR_PROGRAMMING;
			break;
		case VALUE_ITERATION:
			reachabilityMethod = MDPSolnMethod.VALUE_ITERATION;
			break;
		case POLICY_ITERATION:
			reachabilityMethod = MDPSolnMethod.POLICY_ITERATION;
			break;
		case DEFAULT:
		default:
			reachabilityMethod = MDPSolnMethod.GAUSS_SEIDEL;
			break;
		}

		ModelCheckerResult collapsedRes = computeReachProbs(collapsedMdp, null, goal, min, null, null, reachabilityMethod);

		res.timePre += preTime.elapsedSeconds();
		res.timeTaken = time.elapsedSeconds();
		res.numIters += collapsedRes.numIters;
		res.timePre += collapsedRes.timePre;
		res.timeProb0 += collapsedRes.timeProb0;

		// Return results
		for (int i = 0; i < model.getNumStates(); i++) {
			res.soln[i] = maxMECReward * collapsedRes.soln[collapsedMdp.getCollapsedStateIndexForOriginalState(i)];
		}
		return res;
	}

	IntCollection getOptimalChoices(int state, int currentChoice, IntIterator availableChoices, StateChoiceValueFunction function,
			boolean min, boolean all)
	{
		if (!availableChoices.hasNext()) {
			return IntLists.EMPTY_LIST;
		}

		double optimalValue = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
		double currentChoiceValue = Double.NaN;
		Int2DoubleMap choiceValues = all ? new Int2DoubleLinkedOpenHashMap() : null;

		int optimalChoice = -1;
		while (availableChoices.hasNext()) {
			int choice = availableChoices.nextInt();
			double value = function.getValue(state, choice);
			assert !Double.isNaN(value);
			if (all) {
				choiceValues.put(choice, value);
			}
			if (min && value < optimalValue || !min && value > optimalValue) {
				optimalChoice = choice;
				optimalValue = value;
			}
			if (choice == currentChoice) {
				currentChoiceValue = value;
			}
		}

		assert optimalChoice >= 0 && !Double.isNaN(currentChoiceValue) && Double.isFinite(optimalValue);

		if (all) {
			IntList optimalChoices = new IntArrayList();
			for (Int2DoubleMap.Entry entry : choiceValues.int2DoubleEntrySet()) {
				if (PrismUtils.doublesAreCloseAbs(entry.getDoubleValue(), optimalValue, PrismUtils.epsilonDouble)) {
					optimalChoices.add(entry.getIntKey());
				}
			}
			assert optimalChoices.contains(currentChoice) == PrismUtils.doublesAreCloseAbs(currentChoiceValue, optimalValue, PrismUtils.epsilonDouble);
			return optimalChoices;
		}

		if (optimalChoice == currentChoice || PrismUtils.doublesAreCloseAbs(currentChoiceValue, optimalValue, PrismUtils.epsilonDouble)) {
			// The current choice is optimal
			return IntLists.singleton(currentChoice);
		}
		return IntLists.singleton(optimalChoice);
	}

	protected AverageRewardPISolutionMethod getAverageRewardPISolutionMethod() throws PrismException
	{
		String solutionMethodString = settings.getString(PrismSettings.PRISM_MDP_MP_PI_METHOD);
		AverageRewardPISolutionMethod solutionMethod;
		if (solutionMethodString.equals("Full")) {
			solutionMethod = AverageRewardPISolutionMethod.FULL;
		} else if (solutionMethodString.equals("Interleaved bias update")) {
			solutionMethod = AverageRewardPISolutionMethod.FULL_INTERLEAVED_BIAS_UPDATE;
		} else if (solutionMethodString.equals("Attractor")) {
			solutionMethod = AverageRewardPISolutionMethod.ATTRACTOR;
		} else if (solutionMethodString.equals("Gain heuristics")) {
			solutionMethod = AverageRewardPISolutionMethod.GAIN_HEURISTICS;
		} else if (solutionMethodString.equals("Hybrid")) {
			solutionMethod = AverageRewardPISolutionMethod.VALUE_ITERATION_HYBRID;
		} else {
			throw new PrismException("Unknown solution method " + solutionMethodString);
		}
		return solutionMethod;
	}

	/**
	 * Compute expected long run average reward using Policy Iteration
	 */
	public ProcessingResult<? extends GainResult> computeAverageRewardPI(MDP model, MDPRewards mdpRewards, boolean min, MEC mec, MDStrategy policy)
			throws PrismException
	{
		return createPISolver(model, mdpRewards, min, mec, policy, termCritParam, getAverageRewardPISolutionMethod()).solve();
	}

	/**
	 * Compute reachability probabilities using Gauss-Seidel (including Jacobi-style updates).
	 *
	 * @param mdp   The MDP
	 * @param no    Probability 0 states
	 * @param yes   Probability 1 states
	 * @param min   Min or max probabilities (true=min, false=max)
	 * @param init  Optionally, an initial solution vector (will be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 *              Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachProbsGaussSeidel(MDP mdp, NatBitSet no, NatBitSet yes, boolean min, double init[], IntSet known, MDStrategy strat)
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
		mainLog.println("Starting Gauss-Seidel (" + (min ? "min" : "max") + ")...");

		// Store num states
		n = mdp.getNumStates();

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
			maxDiff = mdp.mvMultGSMinMax(soln, min, unknown, termCrit == TermCrit.ABSOLUTE, strat);
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
	 * Compute reachability probabilities using policy iteration.
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param mdp:  The MDP
	 * @param no:   Probability 0 states
	 * @param yes:  Probability 1 states
	 * @param min:  Min or max probabilities (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	protected ModelCheckerResult computeReachProbsPolIter(MDP mdp, NatBitSet no, NatBitSet yes, boolean min, MDStrategy strat) throws PrismException
	{
		ModelCheckerResult res;
		int i, n, iterations, totalIters;
		double soln[], soln2[];
		boolean done;
		long timer;
		DTMCModelChecker mcDTMC;
		DTMC dtmc;

		// Re-use solution to solve each new policy (strategy)?
		boolean reUseSoln = true;

		// Start policy iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting policy iteration (" + (min ? "min" : "max") + ")...");

		// Create a DTMC model checker (for solving policies)
		mcDTMC = new DTMCModelChecker(this);
		mcDTMC.inheritSettings(this);
		mcDTMC.setLog(new PrismDevNullLog());

		// Store num states
		n = mdp.getNumStates();

		// Create solution vectors
		soln = new double[n];
		soln2 = new double[n];

		// Initialise solution vectors.
		for (i = 0; i < n; i++)
			soln[i] = soln2[i] = yes.contains(i) ? 1.0 : 0.0;

		// If not passed in, create new storage for strategy and initialise
		// Initial strategy just picks first choice (0) everywhere
		if (strat == null) {
			strat = new MDStrategyArray(mdp, new int[n]);
			for (i = 0; i < n; i++) {
				strat.setChoiceIndex(i, 0);
			}
		}
		// Otherwise, just initialise for states not in yes/no
		// (Optimal choices for yes/no should already be known)
		else {
			for (i = 0; i < n; i++)
				if (!(no.contains(i) || yes.contains(i)))
					strat.setChoiceIndex(i, 0);
		}

		NatBitSet unknown = NatBitSets.boundedFilledSet(mdp.getNumStates());
		unknown.andNot(no);
		unknown.andNot(yes);

		// Start iterations
		iterations = totalIters = 0;
		done = false;
		while (!done && iterations < maxIters) {
			iterations++;
			// Solve induced DTMC for strategy
			dtmc = new DTMCFromMDPAndMDStrategy(mdp, strat);
			res = mcDTMC.computeReachProbs(dtmc, unknown, yes, reUseSoln ? soln : null, null);
			soln = res.soln;
			totalIters += res.numIters;
			// Check if optimal, improve non-optimal choices
			mdp.mvMultMinMax(soln, min, soln2, NatBitSets.fullSet(mdp.getNumStates()), null);
			done = true;
			for (i = 0; i < n; i++) {
				// Don't look at no/yes states - we may not have strategy info for them,
				// so they might appear non-optimal
				if (no.contains(i) || yes.contains(i))
					continue;
				if (!PrismUtils.doublesAreEqual(soln[i], soln2[i])) {
					done = false;
					IntList opt = mdp.mvMultMinMaxSingleChoices(i, soln, min, soln2[i]);
					// Only update strategy if strictly better
					if (!opt.contains(strat.getChoiceIndex(i))) {
						strat.setChoiceIndex(i, opt.getInt(0));
					}
				}
			}
		}

		// Finished policy iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Policy iteration");
		mainLog.println(" took " + iterations + " cycles (" + totalIters + " iterations in total) and " + timer / 1000.0 + " seconds.");

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iterations + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Return results
		// (Note we don't add the strategy - the one passed in is already there
		// and might have some existing choices stored for other states).
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = totalIters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute reachability probabilities using modified policy iteration.
	 *
	 * @param mdp:  The MDP
	 * @param no:   Probability 0 states
	 * @param yes:  Probability 1 states
	 * @param min:  Min or max probabilities (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	protected ModelCheckerResult computeReachProbsModPolIter(MDP mdp, NatBitSet no, NatBitSet yes, boolean min, MDStrategy strat) throws PrismException
	{
		ModelCheckerResult res;
		int i, n, iters, totalIters;
		double soln[], soln2[];
		boolean done;
		long timer;
		DTMCModelChecker mcDTMC;
		DTMC dtmc;

		// Start value iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting modified policy iteration (" + (min ? "min" : "max") + ")...");

		// Create a DTMC model checker (for solving policies)
		mcDTMC = new DTMCModelChecker(this);
		mcDTMC.inheritSettings(this);
		mcDTMC.setLog(new PrismDevNullLog());

		// Limit iters for DTMC solution - this implements "modified" policy iteration
		mcDTMC.setMaxIters(100);
		mcDTMC.setErrorOnNonConverge(false);

		// Store num states
		n = mdp.getNumStates();

		// Create solution vectors
		soln = new double[n];
		soln2 = new double[n];

		// Initialise solution vectors.
		for (i = 0; i < n; i++)
			soln[i] = soln2[i] = yes.contains(i) ? 1.0 : 0.0;

		// If not passed in, create new storage for strategy and initialise
		// Initial strategy just picks first choice (0) everywhere
		if (strat == null) {
			strat = new MDStrategyArray(mdp, new int[n]);
			for (i = 0; i < n; i++)
				strat.setChoiceIndex(i, 0);
		}
		// Otherwise, just initialise for states not in yes/no
		// (Optimal choices for yes/no should already be known)
		else {
			for (i = 0; i < n; i++)
				if (!(no.contains(i) || yes.contains(i)))
					strat.setChoiceIndex(i, 0);
		}

		// Start iterations
		iters = totalIters = 0;
		done = false;
		while (!done) {
			iters++;
			// Solve induced DTMC for strategy
			dtmc = new DTMCFromMDPAndMDStrategy(mdp, strat);
			res = mcDTMC.computeReachProbsGaussSeidel(dtmc, no, yes, soln, null);
			soln = res.soln;
			totalIters += res.numIters;
			// Check if optimal, improve non-optimal choices
			mdp.mvMultMinMax(soln, min, soln2, NatBitSets.fullSet(mdp.getNumStates()), null);
			done = true;
			for (i = 0; i < n; i++) {
				// Don't look at no/yes states - we don't store strategy info for them,
				// so they might appear non-optimal
				if (no.contains(i) || yes.contains(i))
					continue;
				if (!PrismUtils.doublesAreClose(soln[i], soln2[i], termCritParam, termCrit == TermCrit.ABSOLUTE)) {
					done = false;
					IntList opt = mdp.mvMultMinMaxSingleChoices(i, soln, min, soln2[i]);
					strat.setChoiceIndex(i, opt.getInt(0));
				}
			}
		}

		// Finished policy iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Modified policy iteration");
		mainLog.println(" took " + iters + " cycles (" + totalIters + " iterations in total) and " + timer / 1000.0 + " seconds.");

		// Return results
		// (Note we don't add the strategy - the one passed in is already there
		// and might have some existing choices stored for other states).
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = totalIters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Construct strategy information for min/max reachability probabilities.
	 * (More precisely, list of indices of choices resulting in min/max.)
	 * (Note: indices are guaranteed to be sorted in ascending order.)
	 *
	 * @param mdp      The MDP
	 * @param state    The state to generate strategy info for
	 * @param target   The set of target states to reach
	 * @param min      Min or max probabilities (true=min, false=max)
	 * @param lastSoln Vector of values from which to recompute in one iteration
	 */
	public IntList probReachStrategy(MDP mdp, int state, NatBitSet target, boolean min, double lastSoln[]) throws PrismException
	{
		double val = mdp.mvMultMinMaxSingle(state, lastSoln, min, null);
		return mdp.mvMultMinMaxSingleChoices(state, lastSoln, min, val);
	}

	/**
	 * Compute bounded reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target} within k steps.
	 *
	 * @param mdp    The MDP
	 * @param target Target states
	 * @param k      Bound
	 * @param min    Min or max probabilities (true=min, false=max)
	 */
	public ModelCheckerResult computeBoundedReachProbs(MDP mdp, NatBitSet target, int k, boolean min) throws PrismException
	{
		return computeBoundedReachProbs(mdp, null, target, k, min, null, null);
	}

	/**
	 * Compute bounded until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 *
	 * @param mdp    The MDP
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k      Bound
	 * @param min    Min or max probabilities (true=min, false=max)
	 */
	public ModelCheckerResult computeBoundedUntilProbs(MDP mdp, NatBitSet remain, NatBitSet target, int k, boolean min) throws PrismException
	{
		return computeBoundedReachProbs(mdp, remain, target, k, min, null, null);
	}

	/**
	 * Compute bounded reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 *
	 * @param mdp     The MDP
	 * @param remain  Remain in these states (optional: null means "all")
	 * @param target  Target states
	 * @param k       Bound
	 * @param min     Min or max probabilities (true=min, false=max)
	 * @param init    Optionally, an initial solution vector (may be overwritten)
	 * @param results Optional array of size k+1 to store (init state) results for each step (null if unused)
	 */
	public ModelCheckerResult computeBoundedReachProbs(MDP mdp, NatBitSet remain, NatBitSet target, int k, boolean min, double init[], double results[])
			throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		long timer;

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting bounded probabilistic reachability (" + (min ? "min" : "max") + ")...");

		// Store num states
		n = mdp.getNumStates();

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
			results[0] = Utils.minMaxOverArraySubset(soln2, mdp.getInitialStates(), true);
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
			// Matrix-vector multiply and min/max ops
			mdp.mvMultMinMax(soln, min, soln2, unknown, null);
			// Store intermediate results if required
			// (compute min/max value over initial states for this step)
			if (results != null) {
				// TODO: whether this is min or max should be specified somehow
				results[iters] = Utils.minMaxOverArraySubset(soln2, mdp.getInitialStates(), true);
			}
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Bounded probabilistic reachability (" + (min ? "min" : "max") + ")");
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
	 * Compute expected cumulative (step-bounded) rewards.
	 * i.e. compute the min/max reward accumulated within {@code k} steps.
	 *
	 * @param mdp        The MDP
	 * @param mdpRewards The rewards
	 * @param k          Iteration bound
	 * @param min        Min or max rewards (true=min, false=max)
	 */
	public ModelCheckerResult computeCumulativeRewards(MDP mdp, MDPRewards mdpRewards, int k, boolean min) throws PrismException
	{
		ModelCheckerResult res;
		int i, n, iters;
		long timer;
		double soln[], soln2[], tmpsoln[];

		// Start expected cumulative reward
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting expected cumulative reward (" + (min ? "min" : "max") + ")...");

		// Store num states
		n = mdp.getNumStates();

		// Create/initialise solution vector(s)
		soln = new double[n];
		soln2 = new double[n];
		Arrays.fill(soln, 0.0d);
		Arrays.fill(soln2, 0.0d);

		// Start iterations
		iters = 0;
		while (iters < k) {
			iters++;
			// Matrix-vector multiply and min/max ops
			mdp.mvMultRewMinMax(soln, mdpRewards, min, soln2, NatBitSets.fullSet(mdp.getNumStates()), null);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Expected cumulative reward (" + (min ? "min" : "max") + ")");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;

		return res;
	}

	/**
	 * Compute expected reachability rewards.
	 *
	 * @param mdp        The MDP
	 * @param mdpRewards The rewards
	 * @param target     Target states
	 * @param min        Min or max rewards (true=min, false=max)
	 */
	public ModelCheckerResult computeReachRewards(MDP mdp, MDPRewards mdpRewards, NatBitSet target, boolean min) throws PrismException
	{
		return computeReachRewards(mdp, mdpRewards, target, min, null, null);
	}

	/**
	 * Compute expected reachability rewards.
	 * i.e. compute the min/max reward accumulated to reach a state in {@code target}.
	 *
	 * @param mdp        The MDP
	 * @param mdpRewards The rewards
	 * @param target     Target states
	 * @param min        Min or max rewards (true=min, false=max)
	 * @param init       Optionally, an initial solution vector (may be overwritten)
	 * @param known      Optionally, a set of states for which the exact answer is known
	 *                   Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values).
	 *                   Also, 'known' values cannot be passed for some solution methods, e.g. policy iteration.
	 */
	public ModelCheckerResult computeReachRewards(MDP mdp, MDPRewards mdpRewards, NatBitSet target, boolean min, double init[], NatBitSet known)
			throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet inf;
		int n, numTarget, numInf;
		long timer, timerProb1;
		MDStrategy strat = null;
		// Local copy of setting
		MDPSolnMethod mdpSolnMethod = this.mdpSolnMethod;

		// Switch to a supported method, if necessary
		if (!(mdpSolnMethod == MDPSolnMethod.VALUE_ITERATION || mdpSolnMethod == MDPSolnMethod.GAUSS_SEIDEL
				|| mdpSolnMethod == MDPSolnMethod.POLICY_ITERATION)) {
			mdpSolnMethod = MDPSolnMethod.GAUSS_SEIDEL;
			mainLog.printWarning("Switching to MDP solution method \"" + mdpSolnMethod.fullName() + "\"");
		}

		// Check for some unsupported combinations
		if (mdpSolnMethod == MDPSolnMethod.POLICY_ITERATION) {
			if (known != null) {
				throw new PrismException("Policy iteration methods cannot be passed 'known' values for some states");
			}
		}

		// Start expected reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting expected reachability (" + (min ? "min" : "max") + ")...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		mdp.checkForDeadlocks(target);

		// Store num states
		n = mdp.getNumStates();
		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null && !known.isEmpty()) {
			NatBitSet targetNew = target.clone();
			known.forEach((IntConsumer) i -> {
				if (init[i] == 1.0) {
					targetNew.set(i);
				}
			});
			target = targetNew;
		}

		// If required, export info about target states
		if (getExportTarget()) {
			NatBitSet bsInit = NatBitSets.boundedSet(n);
			for (int i = 0; i < n; i++) {
				bsInit.set(i, mdp.isInitialState(i));
			}
			List<NatBitSet> labels = Arrays.asList(bsInit, target);
			List<String> labelNames = Arrays.asList("init", "target");
			mainLog.println("\nExporting target states info to file \"" + getExportTargetFilename() + "\"...");
			exportLabels(mdp, labels, labelNames, Prism.EXPORT_PLAIN, new PrismFileLog(getExportTargetFilename()));
		}

		// If required, create/initialise strategy storage
		// Set choices to -1, denoting unknown
		// (except for target states, which are -2, denoting arbitrary)
		if (genStrat || exportAdv || mdpSolnMethod == MDPSolnMethod.POLICY_ITERATION) {
			strat = new MDStrategyArray(mdp, new int[n]);
			for (int i = 0; i < n; i++) {
				strat.setChoiceIndex(i, target.contains(i) ? -2 : -1);
			}
		}

		// Precomputation (not optional)
		timerProb1 = System.currentTimeMillis();
		inf = prob1(mdp, null, target, !min, strat);
		inf.flip(0, n);
		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numTarget = target.size();
		numInf = inf.size();
		mainLog.println("target=" + numTarget + ", inf=" + numInf + ", rest=" + (n - (numTarget + numInf)));

		// If required, generate strategy for "inf" states.
		if (genStrat || exportAdv || mdpSolnMethod == MDPSolnMethod.POLICY_ITERATION) {
			assert strat != null;
			if (min) {
				// If min reward is infinite, all choices give infinity
				// So the choice can be arbitrary, denoted by -2;
				IntIterator iterator = inf.iterator();
				while (iterator.hasNext()) {
					strat.setChoiceIndex(iterator.nextInt(), -2);
				}
			} else {
				// If max reward is infinite, there is at least one choice giving infinity.
				// So we pick, for all "inf" states, the first choice for which some transitions stays in "inf".
				IntIterator iterator = inf.iterator();
				while (iterator.hasNext()) {
					int i = iterator.nextInt();
					int numChoices = mdp.getNumChoices(i);
					for (int k = 0; k < numChoices; k++) {
						if (mdp.someSuccessorsInSet(i, k, inf)) {
							strat.setChoiceIndex(i, k);
							break;
						}
					}
				}
			}
		}

		// Compute rewards
		switch (mdpSolnMethod) {
		case VALUE_ITERATION:
			res = computeReachRewardsValIter(mdp, mdpRewards, target, inf, min, init, known, strat);
			break;
		case GAUSS_SEIDEL:
			res = computeReachRewardsGaussSeidel(mdp, mdpRewards, target, inf, min, init, known, strat);
			break;
		case POLICY_ITERATION:
			res = computeReachRewardsPolIter(mdp, mdpRewards, target, inf, min, strat);
			break;
		default:
			throw new PrismException("Unknown MDP solution method " + mdpSolnMethod.fullName());
		}

		// Store strategy
		if (genStrat) {
			res.strat = strat;
		}
		// Export adversary
		if (exportAdv) {
			// Prune strategy
			restrictStrategyToReachableStates(mdp, strat);
			// Export
			PrismLog out = new PrismFileLog(exportAdvFilename);
			new DTMCFromMDPAndMDStrategy(mdp, strat).exportToPrismExplicitTra(out);
			out.close();
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
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param mdp        The MDP
	 * @param mdpRewards The rewards
	 * @param target     Target states
	 * @param inf        States for which reward is infinite
	 * @param min        Min or max rewards (true=min, false=max)
	 * @param init       Optionally, an initial solution vector (will be overwritten)
	 * @param known      Optionally, a set of states for which the exact answer is known
	 * @param strat      Storage for (memoryless) strategy choice indices (ignored if null)
	 *                   Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachRewardsValIter(MDP mdp, MDPRewards mdpRewards, NatBitSet target, NatBitSet inf, boolean min, double init[],
			NatBitSet known,
			MDStrategy strat)
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
		mainLog.println("Starting value iteration (" + (min ? "min" : "max") + ")...");

		// Store num states
		n = mdp.getNumStates();

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
			mdp.mvMultRewMinMax(soln, mdpRewards, min, soln2, unknown, strat);
			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Value iteration (" + (min ? "min" : "max") + ")");
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
	 * Compute expected reachability rewards using Gauss-Seidel (including Jacobi-style updates).
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param mdp        The MDP
	 * @param mdpRewards The rewards
	 * @param target     Target states
	 * @param inf        States for which reward is infinite
	 * @param min        Min or max rewards (true=min, false=max)
	 * @param init       Optionally, an initial solution vector (will be overwritten)
	 * @param known      Optionally, a set of states for which the exact answer is known
	 * @param strat      Storage for (memoryless) strategy choice indices (ignored if null)
	 *                   Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachRewardsGaussSeidel(MDP mdp, MDPRewards mdpRewards, NatBitSet target, NatBitSet inf, boolean min, double init[],
			NatBitSet known, MDStrategy strat) throws PrismException
	{
		ModelCheckerResult res;
		NatBitSet unknown;
		int i, n, iters;
		double soln[], maxDiff;
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting Gauss-Seidel (" + (min ? "min" : "max") + ")...");

		// Store num states
		n = mdp.getNumStates();

		// Create solution vector(s)
		soln = (init == null) ? new double[n] : init;

		// Initialise solution vector. Use (where available) the following in order of preference:
		// (1) exact answer, if already known; (2) 0.0/infinity if in target/inf; (3) passed in initial value; (4) 0.0
		if (init != null) {
			if (known != null) {
				for (i = 0; i < n; i++)
					soln[i] = known.contains(i) ? init[i] : target.contains(i) ? 0.0 : inf.contains(i) ? Double.POSITIVE_INFINITY : init[i];
			} else {
				for (i = 0; i < n; i++)
					soln[i] = target.contains(i) ? 0.0 : inf.contains(i) ? Double.POSITIVE_INFINITY : init[i];
			}
		} else {
			for (i = 0; i < n; i++)
				soln[i] = target.contains(i) ? 0.0 : inf.contains(i) ? Double.POSITIVE_INFINITY : 0.0;
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
			maxDiff = mdp.mvMultRewGSMinMax(soln, mdpRewards, min, unknown, termCrit == TermCrit.ABSOLUTE, strat);
			// Check termination
			done = maxDiff < termCritParam;
		}

		// Finished Gauss-Seidel
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Gauss-Seidel (" + (min ? "min" : "max") + ")");
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
	 * Compute expected reachability rewards using policy iteration.
	 * The array {@code strat} is used both to pass in the initial strategy for policy iteration,
	 * and as storage for the resulting optimal strategy (if needed).
	 * Passing in an initial strategy is required when some states have infinite reward,
	 * to avoid the possibility of policy iteration getting stuck on an infinite-value strategy.
	 *
	 * @param mdp        The MDP
	 * @param mdpRewards The rewards
	 * @param target     Target states
	 * @param inf        States for which reward is infinite
	 * @param min        Min or max rewards (true=min, false=max)
	 * @param strat      Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	protected ModelCheckerResult computeReachRewardsPolIter(MDP mdp, MDPRewards mdpRewards, NatBitSet target, NatBitSet inf, boolean min, MDStrategy strat)
			throws PrismException
	{
		ModelCheckerResult res;
		int i, n, iters, totalIters;
		double soln[], soln2[];
		boolean done;
		long timer;
		DTMCModelChecker mcDTMC;
		DTMC dtmc;
		MCRewards mcRewards;

		// Re-use solution to solve each new policy (strategy)?
		boolean reUseSoln = true;

		// Start policy iteration
		timer = System.currentTimeMillis();
		mainLog.println("Starting policy iteration (" + (min ? "min" : "max") + ")...");

		// Create a DTMC model checker (for solving policies)
		mcDTMC = new DTMCModelChecker(this);
		mcDTMC.inheritSettings(this);
		mcDTMC.setLog(new PrismDevNullLog());

		// Store num states
		n = mdp.getNumStates();

		// Create solution vector(s)
		soln = new double[n];
		soln2 = new double[n];

		// Initialise solution vectors.
		for (i = 0; i < n; i++) {
			double value;
			if (target.contains(i)) {
				value = 0.0;
			} else if (inf.contains(i)) {
				value = Double.POSITIVE_INFINITY;
			} else {
				value = 0.0;
			}
			soln[i] = value;
			soln2[i] = value;
		}

		// If not passed in, create new storage for strategy and initialise
		// Initial strategy just picks first choice (0) everywhere
		if (strat == null) {
			strat = new MDStrategyArray(mdp, new int[n]);
		}

		// Start iterations
		iters = totalIters = 0;
		done = false;
		while (!done && iters < maxIters) {
			iters++;
			// Solve induced DTMC for strategy
			dtmc = new DTMCFromMDPAndMDStrategy(mdp, strat);
			mcRewards = new MCRewardsFromMDPRewards(mdpRewards, strat);
			res = mcDTMC.computeReachRewardsValIter(dtmc, mcRewards, target, inf, reUseSoln ? soln : null, null);
			soln = res.soln;
			totalIters += res.numIters;
			// Check if optimal, improve non-optimal choices
			mdp.mvMultRewMinMax(soln, mdpRewards, min, soln2, NatBitSets.fullSet(mdp.getNumStates()), null);
			done = true;
			for (i = 0; i < n; i++) {
				// Don't look at target/inf states - we may not have strategy info for them,
				// so they might appear non-optimal
				if (target.contains(i) || inf.contains(i)) {
					continue;
				}
				if (!PrismUtils.doublesAreClose(soln[i], soln2[i], termCritParam, termCrit == TermCrit.ABSOLUTE)) {
					done = false;
					IntList opt = mdp.mvMultRewMinMaxSingleChoices(i, soln, mdpRewards, min, soln2[i]);
					// Only update strategy if strictly better
					if (!opt.contains(strat.getChoiceIndex(i))) {
						strat.setChoiceIndex(i, opt.getInt(0));
					}
				}
			}
		}

		// Finished policy iteration
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Policy iteration");
		mainLog.println(" took " + iters + " cycles (" + totalIters + " iterations in total) and " + timer / 1000.0 + " seconds.");

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = totalIters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Construct strategy information for min/max expected reachability.
	 * (More precisely, list of indices of choices resulting in min/max.)
	 * (Note: indices are guaranteed to be sorted in ascending order.)
	 *
	 * @param mdp        The MDP
	 * @param mdpRewards The rewards
	 * @param state      The state to generate strategy info for
	 * @param target     The set of target states to reach
	 * @param min        Min or max rewards (true=min, false=max)
	 * @param lastSoln   Vector of values from which to recompute in one iteration
	 */
	public IntList expReachStrategy(MDP mdp, MDPRewards mdpRewards, int state, NatBitSet target, boolean min, double lastSoln[]) throws PrismException
	{
		double val = mdp.mvMultRewMinMaxSingle(state, lastSoln, mdpRewards, min, null);
		return mdp.mvMultRewMinMaxSingleChoices(state, lastSoln, mdpRewards, min, val);
	}

	/**
	 * Restrict a (memoryless) strategy for an MDP, stored as an integer array of choice indices,
	 * to the states of the MDP that are reachable under that strategy.
	 *
	 * @param mdp   The MDP
	 * @param strat The strategy
	 */
	public void restrictStrategyToReachableStates(MDP mdp, MDStrategy strat)
	{
		NatBitSet restrict = NatBitSets.boundedSet(mdp.getNumStates());
		NatBitSet explore = NatBitSets.boundedSet(mdp.getNumStates());
		// Get initial states
		mdp.getInitialStates().forEach((IntConsumer) initialState -> {
			restrict.set(initialState);
			explore.set(initialState);
		});
		// Compute reachable states (store in 'restrict')
		boolean foundMore = true;
		while (foundMore) {
			foundMore = false;
			IntIterator iterator = explore.iterator();
			while (iterator.hasNext()) {
				int s = iterator.nextInt();
				explore.set(s, false);
				if (strat.isChoiceDefined(s)) {
					for (Int2DoubleMap.Entry e : mdp.getTransitions(s, strat.getChoiceIndex(s))) {
						int dest = e.getIntKey();
						if (!restrict.contains(dest)) {
							foundMore = true;
							restrict.set(dest);
							explore.set(dest);
						}
					}
				}
			}
		}
		// Set strategy choice for non-reachable state to -1
		int n = mdp.getNumStates();
		NatBitSets.complementIterator(restrict, n).forEachRemaining((IntConsumer) s -> strat.setChoiceIndex(s, -3));
	}

	private PolicyIterationSolver createPISolver(MDP model, MDPRewards mdpRewards, boolean min, MEC mec, MDStrategy policy, double requiredPrecision,
			AverageRewardPISolutionMethod solutionMethod) throws PrismException
	{
		if ((solutionMethod == AverageRewardPISolutionMethod.ATTRACTOR || solutionMethod == AverageRewardPISolutionMethod.VALUE_ITERATION_HYBRID)
				&& mec == null) {
			if (model.getMDPType() == MDP.Type.UNKNOWN) {
				mainLog.printWarning("Attractor PI and modified PI only work on MECs - assuming that the model is communicating");
			} else if (model.getMDPType() == MDP.Type.MULTICHAIN) {
				throw new PrismException("Attractor PI and modified PI do not work for multi chain models");
			}
		}

		if (solutionMethod == AverageRewardPISolutionMethod.VALUE_ITERATION_HYBRID) {
			return new PolicyIterationHybridSolver(model, mdpRewards, min, mec, policy, requiredPrecision, 10);
		}
		if (solutionMethod == AverageRewardPISolutionMethod.GAIN_HEURISTICS) {
			return new PolicyIterationHeuristicSolver(model, mdpRewards, min, mec, policy, solutionMethod, requiredPrecision);
		}
		return new PolicyIterationFullSolver(model, mdpRewards, min, mec, policy, solutionMethod, AverageRewardPIEvaluationMethod.LINEAR_EQUATION);
	}

	public enum AverageRewardPIMECEvaluationMethod
	{
		FULL, SCC_DECOMPOSITION
	}

	public enum AverageRewardPIEvaluationMethod
	{
		LINEAR_EQUATION, VALUE_ITERATION_HEURISTIC, VALUE_ITERATION
	}

	public enum AverageRewardPISolutionMethod
	{
		FULL, FULL_INTERLEAVED_BIAS_UPDATE, VALUE_ITERATION_HYBRID, GAIN_HEURISTICS, ATTRACTOR
	}

	public enum AverageRewardMECSolutionMethod
	{
		VALUE_ITERATION, POLICY_ITERATION, MODIFIED_POLICY_ITERATION, ATTRACTOR_POLICY_ITERATION, INTERVAL_POLICY_ITERATION
	}

	public enum AverageRewardMECReachMethod
	{
		INTERVAL_VALUE_ITERATION, LINEAR_PROGRAM, VALUE_ITERATION, POLICY_ITERATION, DEFAULT
	}

	public interface MecAverageRewardSolverFactory
	{
		MecAverageRewardSolver obtain(MDP mdp, MDPRewards rewards, PrismSettings settings, boolean min, double precision);
	}

	@FunctionalInterface
	public interface MecAverageRewardSolver
	{
		MeanPayoffModelCheckerResult solve(MEC mec) throws PrismException;
	}

	@FunctionalInterface
	private interface StateChoiceValueFunction
	{
		double getValue(int state, int choice);
	}

	public static class MeanPayoffModelCheckerResult extends ModelCheckerResult
	{
		public double meanPayoff;

		@Override public void clear()
		{
			super.clear();
			meanPayoff = 0;
		}
	}

	public static class ValueIterationModelCheckerResult extends MeanPayoffModelCheckerResult
	{
		public double maximumValue;
		public double minimumValue;

		@Override public void clear()
		{
			super.clear();
			maximumValue = 0d;
			minimumValue = 0d;
		}
	}

	public static class MDPSpanNormValueIterator extends SpanNormValueIterator
	{
		private final MDPExplicit model;
		private final MDPRewards rewards;
		private final boolean min;

		MDPSpanNormValueIterator(PrismComponent parent, MDPExplicit model, MDPRewards rewards, boolean min)
		{
			super(parent);
			this.model = model;
			this.rewards = rewards;
			this.min = min;
		}

		@Override MDPExplicit getModel()
		{
			return model;
		}

		public MDPRunner createRunner(MEC mec, double precision)
		{
			return new MDPRunner(mec, step -> (step - 1) % 5 == 0, precision);
		}

		class MDPRunner extends Runner
		{
			private final MEC mec;
			private double resultMeanPayoff;

			MDPRunner(MEC mec, IntPredicate checkStepPredicate, double precision)
			{
				super(checkStepPredicate, precision);
				this.mec = mec;

				String minMaxString = min ? "min" : "max";
				if (mec == null) {
					mainLog.println(String.format("Starting %s value iteration on MDP with %d states", minMaxString, model.getStatesList().size()),
							PrismLog.VL_HIGH);
				} else {
					mainLog.println(String.format("Starting %s value iteration on MDP subset of size %d", minMaxString, mec.states.size()), PrismLog.VL_HIGH);
				}
			}

			@Override double operator(int state, IntToDoubleFunction stateValues, double tau)
			{
				IntSet allowedActions = mec == null ? null : mec.actions.get(state);
				return model.mvMultRewMinMaxAperiodicSingle(state, allowedActions, stateValues, rewards, min, null, tau);
			}

			@Override IntIterator getStateIterator()
			{
				return FastUtils.iterator(mec, model.getNumStates());
			}
		}
	}

	public class PolicyIterationSolverFactory implements MecAverageRewardSolverFactory
	{
		@Override public MecAverageRewardSolver obtain(MDP mdp, MDPRewards rewards, PrismSettings settings, boolean min, double precision)
		{
			return mec -> {
				PolicyIterationSolver solver = createPISolver(mdp, rewards, min, mec, null, precision, getAverageRewardPISolutionMethod());
				ProcessingResult<? extends GainResult> processingResult = solver.solve();
				MeanPayoffModelCheckerResult result = new MeanPayoffModelCheckerResult();
				result.meanPayoff = processingResult.getResult().getGain(mec.states.firstInt());
				result.numIters = processingResult.getNumIters();
				result.timeTaken = processingResult.getTimeTaken();
				return result;
			};
		}
	}

	public class ValueIterationSolverFactory implements MecAverageRewardSolverFactory
	{
		@Override public MecAverageRewardSolver obtain(MDP mdp, MDPRewards rewards, PrismSettings settings, boolean min, double precision)
		{
			boolean enableAdaptiveTau = settings.getBoolean(PrismSettings.PRISM_MDP_MP_ADAPTIVE_TAU);
			double tauValue = settings.getDouble(PrismSettings.PRISM_MDP_MP_APERIODICITYTAU);

			MDPSpanNormValueIterator valueIterator =
					new MDPSpanNormValueIterator(MDPModelChecker.this, (MDPExplicit) mdp, rewards, min);
			return mec -> {
				MDPSpanNormValueIterator.MDPRunner runner = valueIterator.createRunner(mec, precision);
				runner.run(-1);
				MeanPayoffModelCheckerResult result = new MeanPayoffModelCheckerResult();
				result.meanPayoff = runner.getResultMeanPayoff();
				result.numIters = runner.getIterations();
				result.timeTaken = runner.getTimeTaken();
				result.meanPayoff = runner.getResultMeanPayoff();
				return result;
			};
		}
	}

	private class PolicyIterationHybridSolver extends PolicyIterationSolver
	{
		private final int stepsPerPolicy;

		PolicyIterationHybridSolver(MDP model, MDPRewards mdpRewards, boolean min, MEC mec, MDStrategy policy, double requiredPrecision, int stepsPerPolicy)
				throws PrismException
		{
			super(model, mdpRewards, min, mec, policy, requiredPrecision, AverageRewardPISolutionMethod.VALUE_ITERATION_HYBRID);
			this.stepsPerPolicy = stepsPerPolicy;
		}

		@Override ProcessingResult<? extends GainResult> solve()
		{
			Time time = new Time();

			DTMCModelChecker.MCSpanNormValueIterator iterator =
					new DTMCModelChecker.MCSpanNormValueIterator(MDPModelChecker.this, dtmc, mcRewards);
			DTMCModelChecker.MCSpanNormValueIterator.MCRunner runner =
					iterator.createRunner(getSubset(), step -> step % stepsPerPolicy == 0, requiredPrecision);

			StateChoiceValueFunction successorGainFunction = (state, choice) -> model.mvMultSingle(state, choice,
					s -> mdpRewards.getTransitionReward(state, choice) + runner.getValue(s));

			boolean converged = false;
			while (!converged) {
				steps += 1;
				converged = runner.run(stepsPerPolicy);

				// Improve according to total accumulated value
				IntIterator improvementStates = getImprovementStates();
				while (improvementStates.hasNext()) {
					int state = improvementStates.nextInt();
					assert policy.isChoiceDefined(state);
					int currentChoice = policy.getChoiceIndex(state);

					IntIterator allowedChoices = getAllowedActions(state);
					IntCollection optimalChoices = getOptimalChoices(state, currentChoice, allowedChoices, successorGainFunction, min, false);
					policy.setChoiceIndex(state, optimalChoices.iterator().nextInt());
				}
			}

			// Finished policy iteration
			double elapsedSeconds = time.elapsedSeconds();
			mainLog.println(String.format("Modified policy iteration took %d steps, %d iterations, %d improvements, and %g seconds.",
					steps, runner.getIterations(), stateImprovements, elapsedSeconds), PrismLog.VL_HIGH);
			double gain = runner.getResultMeanPayoff();
			assert gain >= 0;
			return new ProcessingResult<>(new GainResultSingle(gain), runner.getObtainedPrecision(), steps + runner.getIterations(), elapsedSeconds);
		}
	}

	private class PolicyIterationHeuristicSolver extends PolicyIterationSolver
	{
		PolicyIterationHeuristicSolver(MDP model, MDPRewards mdpRewards, boolean min, MEC mec, MDStrategy policy,
				AverageRewardPISolutionMethod solutionMethod, double precision) throws PrismException
		{
			super(model, mdpRewards, min, mec, policy, precision, solutionMethod);
		}

		@Override ProcessingResult<? extends GainResult> solve() throws PrismException
		{
			double approximationPrecision = 1e-2;
			int gainImprovements = 0;
			int biasImprovements = 0;

			if (!getImprovementStates().hasNext()) {
				ProcessingResult<IntervalResult> processingResult =
						dtmcModelChecker.computeAverageRewardIntervalValueIteration(dtmc, mcRewards, mec == null ? null : mec.states, requiredPrecision);
				IntervalResult intervalResult = processingResult.getResult();
				GainResult gainResult = new GainResultSparse(intervalResult::getAverage);
				return new ProcessingResult<>(gainResult, processingResult.getPrecision(), steps, time.elapsedSeconds());
			}

			boolean inducedChainIsUnichain = isInducedChainAlwaysUnichain();
			NatBitSet improvedStates = NatBitSets.boundedSet(numMdpStates);
			NatBitSet unknownStates = NatBitSets.boundedSet(numMdpStates);

			Int2ObjectMap<IntCollection> gainOptimalChoices = new Int2ObjectOpenHashMap<>();
			GainResult currentGainResult = null;

			while (true) {
				// We updated the underlying strategy
				dtmc.clearPredecessorRelation();
				gainOptimalChoices.clear();
				steps++;

				// Evaluate policy using value iteration
				// TODO do this with some kind of iteration bound maybe?
				ProcessingResult<IntervalResult> processingResult =
						dtmcModelChecker.computeAverageRewardIntervalValueIteration(dtmc, mcRewards, getSubset(), approximationPrecision);
				if (approximationPrecision > requiredPrecision) {
					approximationPrecision *= 0.9d;
				}
				IntervalResult result = processingResult.getResult();

				// Consistency check
				GainResult previousGainResult = currentGainResult;
				assert previousGainResult == null || Iterators.all(FastUtils.iterator(mec, model.getNumStates()),
						s -> PrismUtils.doublesAreLessOrEqual(previousGainResult.getGain(s), result.getUpper(s)));

				if (mec == null) {
					// Local improvement
					StateChoiceValueFunction successorGainLowerBounds = (state, choice) -> model.mvMultSingle(state, choice, result::getLower);
					StateChoiceValueFunction successorGainUpperBounds = (state, choice) -> model.mvMultSingle(state, choice, result::getUpper);

					doGainImprovementsWithBounds(successorGainLowerBounds, successorGainUpperBounds, gainOptimalChoices, improvedStates, unknownStates);
					boolean improvement = !improvedStates.isEmpty();
					if (improvement) {
						gainImprovements++;
						mainLog.println(String.format("Improving strategy in %d states based on approximate gain estimate", improvedStates.size()),
								PrismLog.VL_HIGH);
					}
					improvedStates.clear();
					unknownStates.clear();
					if (improvement) {
						// We changed something based on the gain evaluation - re-evaluate
						continue;
					}
				} else {
					// Do attractor based improvement
					double maximalLowerBound = 0.0d;
					IntIterator iterator = getSubset().iterator();
					while (iterator.hasNext()) {
						int state = iterator.nextInt();
						double lowerBound = result.getLower(state);
						assert lowerBound <= result.getUpper(state);
						if (lowerBound > maximalLowerBound) {
							maximalLowerBound = lowerBound;
						}
					}
					NatBitSet attractor = NatBitSets.boundedSet(numMdpStates);
					iterator = getSubset().iterator();
					while (iterator.hasNext()) {
						int state = iterator.nextInt();
						double upperBound = result.getUpper(state);
						if (upperBound >= maximalLowerBound) {
							attractor.set(state);
						} else {
							improvedStates.set(state);
						}
					}

					assert !attractor.isEmpty();
					if (!improvedStates.isEmpty()) {
						int improvements = updatePolicyToAttractor(attractor, getModelPredecessorRelation());
						assert 0 < improvements && improvements <= improvedStates.size();
						improvedStates.clear();
						mainLog.println(String.format("Improving strategy in %d states based on approximate gain estimate", improvements), PrismLog.VL_HIGH);
						gainImprovements++;
						continue;
					}
				}

				// Now we evaluate gain and bias precisely
				ProcessingResult<GainBiasResult> gainBiasResult = dtmcModelChecker.computeGainBias(dtmc, mcRewards, getSubset());
				GainBiasResult gainBias = gainBiasResult.getResult();
				currentGainResult = gainBias;

				assert Iterators.all(FastUtils.iterator(mec, model.getNumStates()),
						s -> PrismUtils.doublesAreLessOrEqual(result.getLower(s), gainBias.getGain(s))
								&& PrismUtils.doublesAreLessOrEqual(gainBias.getGain(s), result.getUpper(s)));

				StateChoiceValueFunction successorGainFunction = (state, choice) -> model.mvMultSingle(state, choice, gainBias::getGain);
				int improvements = doGainImprovements(successorGainFunction, gainOptimalChoices);
				if (improvements > 0) {
					mainLog.println(String.format("Improving strategy in %d states based on precise gain", improvements), PrismLog.VL_HIGH);
					gainImprovements++;
					continue;
				}

				StateChoiceValueFunction successorBiasFunction =
						(state, choice) -> model.mvMultSingle(state, choice, gainBias::getBias) + mdpRewards.getTransitionReward(state, choice);
				improvements = doBiasImprovements(successorBiasFunction, gainOptimalChoices);
				if (improvements > 0) {
					mainLog.print(String.format("Improving strategy in %d states based on precise bias", improvements), PrismLog.VL_HIGH);
					biasImprovements++;
					continue;
				}
				break;
			}

			double elapsedSeconds = time.elapsedSeconds();
			mainLog.println(String.format("Heuristic policy iteration took %d steps, %d/%d gain/bias improvement steps, %d improvements, and %g seconds.",
					steps, gainImprovements, biasImprovements, stateImprovements, elapsedSeconds), PrismLog.VL_HIGH);
			return new ProcessingResult<>(currentGainResult, requiredPrecision, steps, elapsedSeconds);
		}
	}

	private class PolicyIterationFullSolver extends PolicyIterationSolver
	{
		private final AverageRewardPIEvaluationMethod evaluationMethod;

		PolicyIterationFullSolver(MDP model, MDPRewards mdpRewards, boolean min, MEC mec, MDStrategy policy,
				AverageRewardPISolutionMethod solutionMethod, AverageRewardPIEvaluationMethod evaluationMethod) throws PrismException
		{
			super(model, mdpRewards, min, mec, policy, PrismUtils.epsilonDouble, solutionMethod);
			this.evaluationMethod = evaluationMethod;
		}

		@Override ProcessingResult<? extends GainResult> solve() throws PrismException
		{
			int gainImprovements = 0;
			int biasImprovements = 0;
			ProcessingResult<GainBiasResult> optimalSolution;
			Int2ObjectMap<IntCollection> gainOptimalChoices = new Int2ObjectOpenHashMap<>(numStates);

			boolean inducedChainIsUnichain = isInducedChainAlwaysUnichain();

			while (true) {
				dtmc.clearPredecessorRelation();
				steps++;

				// Evaluate policy
				ProcessingResult<GainBiasResult> currentSolution;
				if (inducedChainIsUnichain) {
					currentSolution = dtmcModelChecker.computeGainBiasLinearEquationsUnichain(dtmc, mcRewards, getSubset(), null);
				} else {
					currentSolution = dtmcModelChecker.computeGainBias(dtmc, mcRewards, getSubset());
				}

				// Improve based on the result
				GainBiasResult gainBias = currentSolution.getResult();
				boolean gainImprove = false;
				boolean biasImprove = false;

				StateChoiceValueFunction successorGainFunction = (state, choice) -> model.mvMultSingle(state, choice, gainBias::getGain);
				StateChoiceValueFunction successorBiasFunction =
						(state, choice) -> model.mvMultSingle(state, choice, gainBias::getBias) + mdpRewards.getTransitionReward(state, choice);

				if (solutionMethod == AverageRewardPISolutionMethod.FULL || solutionMethod == AverageRewardPISolutionMethod.FULL_INTERLEAVED_BIAS_UPDATE) {
					int improvements = doGainImprovements(successorGainFunction, gainOptimalChoices);
					if (improvements > 0) {
						gainImprovements++;
						gainImprove = true;
					}
					if (!gainImprove || (solutionMethod == AverageRewardPISolutionMethod.FULL_INTERLEAVED_BIAS_UPDATE && !gainOptimalChoices.isEmpty())) {
						// Bias improvement

						// If a state is not contained in gainOptimalChoices, its either a deadlock state or
						// we have interleavedBiasImprovement where this state actually improved his gain
						improvements = doBiasImprovements(successorBiasFunction, gainOptimalChoices);
						if (improvements > 0) {
							biasImprovements++;
							biasImprove = true;
						}
					}
					if (!gainImprove && !biasImprove) {
						// Done
						optimalSolution = currentSolution;
						break;
					}
				} else if (solutionMethod == AverageRewardPISolutionMethod.ATTRACTOR) {
					if (!inducedChainIsUnichain) {
						// Improve gain by attractor
						List<NatBitSet> bsccs = ECComputerFast.computeBSCCs(dtmc, getSubset());
						assert !bsccs.isEmpty();
						// If size == 1, the DTMC is already unichain
						if (bsccs.size() > 1) {
							NatBitSet optimalBscc = null;
							double optimalGain = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
							for (NatBitSet bscc : bsccs) {
								int stateInBscc = bscc.firstInt();
								double bsccGain = gainBias.getGain(stateInBscc);
								if (min && bsccGain < optimalGain || !min && bsccGain > optimalGain) {
									optimalGain = bsccGain;
									optimalBscc = bscc;
								}
							}
							assert optimalBscc != null;
							int improvements = updatePolicyToAttractor(optimalBscc, getModelPredecessorRelation());
							if (improvements > 0) {
								gainImprovements += improvements;
								gainImprove = true;
							}
						}
						inducedChainIsUnichain = true;
					}
					if (!gainImprove) {
						assert inducedChainIsUnichain;
						// Improve bias
						IntIterator stateIterator = getImprovementStates();
						while (stateIterator.hasNext()) {
							int state = stateIterator.nextInt();
							assert policy.isChoiceDefined(state);
							int currentChoice = policy.getChoiceIndex(state);
							// Since this is a unichain model, all states have the same gain, thus every action is allowed
							IntIterator allowedChoices = getAllowedActions(state);
							int optimalBiasChoice =
									getOptimalChoices(state, currentChoice, allowedChoices, successorBiasFunction, min, false).iterator().nextInt();
							if (optimalBiasChoice == currentChoice) {
								continue;
							}
							assert optimalBiasChoice >= 0;
							policy.setChoiceIndex(state, optimalBiasChoice);
							biasImprovements++;
							biasImprove = true;
						}

						if (!biasImprove) {
							optimalSolution = currentSolution;
							break;
						}
						inducedChainIsUnichain = isInducedChainAlwaysUnichain();
					}
				}
			}
			// Finished policy iteration
			double elapsedSeconds = time.elapsedSeconds();
			mainLog.println(String.format("Policy iteration took %d steps, %d/%d gain/bias improvement steps, %d improvements, and %g seconds.",
					steps, gainImprovements, biasImprovements, stateImprovements, elapsedSeconds), PrismLog.VL_HIGH);
			return new ProcessingResult<GainResult>(optimalSolution.getResult(), PrismUtils.epsilonDouble, steps, elapsedSeconds);
		}
	}

	private abstract class PolicyIterationSolver
	{
		final MDP model;
		final MDPRewards mdpRewards;
		final boolean min;
		final MEC mec;
		final MDStrategy policy;
		final AverageRewardPISolutionMethod solutionMethod;
		final DTMCModelChecker dtmcModelChecker;
		final DTMCExplicit dtmc;
		final MCRewards mcRewards;
		final double requiredPrecision;
		final int numMdpStates;
		final int numStates;
		final Time time;
		Supplier<IntIterator> improvementStateSupplier = null;
		int steps = 0;
		int stateImprovements = 0;
		private PredecessorRelation predecessorRelation;

		PredecessorRelation getModelPredecessorRelation() {
			if (predecessorRelation == null) {
				if (mec == null || numMdpStates == numStates) {
					predecessorRelation = model.getPredecessorRelation(MDPModelChecker.this, true);
				} else {
					predecessorRelation = new PredecessorRelationSparse(model, mec.states);
				}
			}
			return predecessorRelation;
		}

		PolicyIterationSolver(MDP model, MDPRewards mdpRewards, boolean min, MEC mec, MDStrategy policy, double requiredPrecision,
				AverageRewardPISolutionMethod solutionMethod) throws PrismException
		{
			this.model = model;
			this.mdpRewards = mdpRewards;
			this.min = min;
			this.mec = mec;
			this.solutionMethod = solutionMethod;

			numMdpStates = model.getNumStates();
			numStates = mec == null ? numMdpStates : mec.size();

			// If not passed, create new storage for strategy and initialise
			// Initial strategy just picks first available choice everywhere
			if (policy == null) {
				if (numMdpStates == numStates) {
					int[] choices = new int[model.getNumStates()];
					this.policy = new MDStrategyArray(model, choices);
				} else {
					Int2IntMap choices = new Int2IntOpenHashMap();
					for (Int2ObjectMap.Entry<NatBitSet> entry : mec.actions.int2ObjectEntrySet()) {
						int state = entry.getIntKey();
						int choice = entry.getValue().firstInt();
						assert choice >= 0;
						choices.put(state, choice);
					}
					choices.defaultReturnValue(-1);
					this.policy = new MDStrategySparse(model, choices);
				}
			} else {
				this.policy = policy;
			}
			if (mec != null) {
				assert Iterators.all(mec.states.iterator(), this.policy::isChoiceDefined);
			}

			// Create a DTMC model checker (for evaluating policies)
			dtmcModelChecker = new DTMCModelChecker(MDPModelChecker.this);
			dtmcModelChecker.inheritSettings(MDPModelChecker.this);
			dtmc = new DTMCFromMDPAndMDStrategy(model, this.policy);
			mcRewards = new MCRewardsFromMDPRewards(mdpRewards, this.policy);
			this.requiredPrecision = requiredPrecision;
			time = new Time();
		}

		IntIterator getImprovementStates()
		{
			if (improvementStateSupplier == null) {
				if (mec == null) {
					NatBitSet deterministicStates = NatBitSets.boundedSet(numMdpStates);
					for (int state = 0; state < model.getNumStates(); state++) {
						if (model.getNumChoices(state) <= 1) {
							deterministicStates.set(state);
						}
					}
					improvementStateSupplier = () -> NatBitSets.complementIterator(deterministicStates, model.getNumStates());
				} else {
					NatBitSet improvementStates = mec.states.clone();
					IntIterator stateIterator = mec.states.iterator();
					while (stateIterator.hasNext()) {
						int state = stateIterator.nextInt();
						if (mec.actions.get(state).size() <= 1) {
							improvementStates.clear(state);
						}
					}
					improvementStateSupplier = improvementStates::iterator;
				}
			}
			return improvementStateSupplier.get();
		}

		boolean isInducedChainAlwaysUnichain()
		{
			return EnumSet.of(MDP.Type.UNICHAIN, MDP.Type.RECURRENT).contains(model.getMDPType());
		}

		boolean isStateConsidered(int state)
		{
			return mec == null || mec.states.contains(state);
		}

		NatBitSet getSubset()
		{
			return mec == null ? null : mec.states;
		}

		IntIterator getAllowedActions(int state)
		{
			if (mec == null) {
				assert model.getNumChoices(state) > 1;
			} else {
				assert mec.actions.get(state).size() > 1;
			}
			return mec == null ? IntIterators.fromTo(0, model.getNumChoices(state)) : mec.actions.get(state).iterator();
		}

		abstract ProcessingResult<? extends GainResult> solve() throws PrismException;

		int doBiasImprovements(StateChoiceValueFunction successorBiasFunction, Int2ObjectMap<IntCollection> stateAllowedChoices)
		{
			int improvements = 0;
			for (Int2ObjectMap.Entry<IntCollection> entry : stateAllowedChoices.int2ObjectEntrySet()) {
				int state = entry.getIntKey();
				assert policy.isChoiceDefined(state);
				int currentChoice = policy.getChoiceIndex(state);
				IntCollection allowedChoices = entry.getValue();

				int optimalBiasChoice =
						getOptimalChoices(state, currentChoice, allowedChoices.iterator(), successorBiasFunction, min, false).iterator().nextInt();
				if (optimalBiasChoice == currentChoice) {
					continue;
				}
				assert optimalBiasChoice >= 0;
				policy.setChoiceIndex(state, optimalBiasChoice);
				improvements++;
			}
			this.stateImprovements += improvements;
			return improvements;
		}

		int doGainImprovements(StateChoiceValueFunction successorGainFunction, Int2ObjectFunction<IntCollection> gainOptimalChoices)
		{
			int improvements = 0;
			IntIterator gainStateIterator = getImprovementStates();
			while (gainStateIterator.hasNext()) {
				int state = gainStateIterator.nextInt();
				assert policy.isChoiceDefined(state);
				int currentChoice = policy.getChoiceIndex(state);
				IntIterator allowedChoices = getAllowedActions(state);
				IntCollection optimalChoices = getOptimalChoices(state, currentChoice, allowedChoices, successorGainFunction, min, true);

				if (optimalChoices.isEmpty()) {
					continue;
				}

				gainOptimalChoices.put(state, optimalChoices);
				if (optimalChoices.contains(currentChoice)) {
					// Did not improve gain
					continue;
				}
				int optimalGainChoice = optimalChoices.iterator().nextInt();
				assert optimalGainChoice >= 0;
				policy.setChoiceIndex(state, optimalGainChoice);
				improvements++;
			}
			this.stateImprovements += improvements;
			return improvements;
		}

		int updatePolicyToAttractor(NatBitSet attractor, PredecessorRelation relation)
		{
			assert getSubset() == null && PredecessorRelations.calculatePreStar(relation, null, attractor, null).size() == numMdpStates
					|| PredecessorRelations.calculatePreStar(relation, getSubset(), attractor, null).equals(getSubset());

			// TODO This can be done faster I think. This essentially is prob1 + strategy generation ...
			int improvements = 0;

			attractor = NatBitSets.ensureModifiable(attractor);
			NatBitSet previouslyAddedStates = NatBitSets.modifiableCopyOf(attractor);
			NatBitSet addedStates = NatBitSets.boundedSet(numMdpStates);
			while (!previouslyAddedStates.isEmpty()) {
				IntIterator iterator = previouslyAddedStates.iterator();
				while (iterator.hasNext()) {
					IntIterator preIterator = relation.getPredecessorsIterator(iterator.nextInt());
					while (preIterator.hasNext()) {
						int predecessor = preIterator.nextInt();
						if (attractor.contains(predecessor)) {
							continue;
						}
						if (!isStateConsidered(predecessor)) {
							continue;
						}
						assert policy.isChoiceDefined(predecessor);
						if (!model.allSuccessorsInSet(predecessor, policy.getChoiceIndex(predecessor), attractor)) {
							int choices = model.getNumChoices(predecessor);

							int attractorChoice = -1;
							for (int choice = 0; choice < choices && attractorChoice == -1; choice++) {
								for (Int2DoubleMap.Entry entry : model.getTransitions(predecessor, choice)) {
									if (attractor.contains(entry.getIntKey())) {
										attractorChoice = choice;
										break;
									}
								}
							}
							assert attractorChoice >= 0;
							policy.setChoiceIndex(predecessor, attractorChoice);
							improvements++;
						}
						attractor.set(predecessor);
						addedStates.set(predecessor);
					}
				}
				assert !previouslyAddedStates.intersects(addedStates);

				NatBitSet swap = previouslyAddedStates;
				previouslyAddedStates = addedStates;
				addedStates = swap;
				addedStates.clear();
			}
			assert attractor.size() == numStates;
			this.stateImprovements += improvements;
			return improvements;
		}

		void doGainImprovementsWithBounds(StateChoiceValueFunction successorLowerBounds, StateChoiceValueFunction successorUpperBounds,
				Int2ObjectMap<IntCollection> optimalActions, NatBitSet improvedStates, NatBitSet unknownStates)
		{
			int improvements = 0;
			StateChoiceValueFunction worstBounds = min ? successorUpperBounds : successorLowerBounds;
			StateChoiceValueFunction bestBounds = min ? successorLowerBounds : successorUpperBounds;
			DoubleComparator comparator = min ? DoubleComparators.OPPOSITE_COMPARATOR : DoubleComparators.NATURAL_COMPARATOR;

			IntIterator gainStateIterator = getImprovementStates();
			Int2DoubleMap choiceBestBounds = new Int2DoubleOpenHashMap();
			while (gainStateIterator.hasNext()) {
				int state = gainStateIterator.nextInt();
				assert policy.isChoiceDefined(state);
				int currentChoice = policy.getChoiceIndex(state);
				IntIterator allowedChoices = getAllowedActions(state);
				assert allowedChoices.hasNext();

				double otherOptimalWorstBound = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
				double otherOptimalBestBound = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
				int otherOptimalWorstBoundChoice = -1;
				int otherOptimalBestBoundChoice = -1;
				double currentChoiceWorstBound = Double.NaN;
				double currentChoiceBestBound = Double.NaN;

				choiceBestBounds.clear();
				while (allowedChoices.hasNext()) {
					int choice = allowedChoices.nextInt();
					double worstBound = worstBounds.getValue(state, choice);
					double bestBound = bestBounds.getValue(state, choice);
					assert comparator.compare(worstBound, bestBound) <= 0;
					if (choice == currentChoice) {
						assert Double.isNaN(currentChoiceWorstBound) && Double.isNaN(currentChoiceBestBound);
						currentChoiceWorstBound = worstBound;
						currentChoiceBestBound = bestBound;
					} else {
						choiceBestBounds.put(state, bestBound);

						if (comparator.compare(worstBound, otherOptimalWorstBound) > 0) {
							assert min && worstBound < otherOptimalWorstBound || !min && worstBound > otherOptimalWorstBound;
							otherOptimalWorstBound = worstBound;
							otherOptimalWorstBoundChoice = choice;
						}
						if (comparator.compare(bestBound, otherOptimalBestBound) > 0) {
							assert min && bestBound < otherOptimalBestBound || !min && bestBound > otherOptimalBestBound;
							otherOptimalBestBound = bestBound;
							otherOptimalBestBoundChoice = choice;
						}
					}
				}
				assert Double.isFinite(otherOptimalBestBound) && Double.isFinite(otherOptimalWorstBound);
				assert otherOptimalBestBoundChoice >= 0 && otherOptimalWorstBoundChoice >= 0;
				assert min && otherOptimalBestBound <= otherOptimalWorstBound ||
						!min && otherOptimalWorstBound <= otherOptimalBestBound;

				NatBitSet stateOptimalActions = NatBitSets.boundedSet(model.getNumChoices(state));
				optimalActions.put(state, stateOptimalActions);
				if (comparator.compare(currentChoiceWorstBound, otherOptimalBestBound) > 0) {
					stateOptimalActions.set(currentChoice);
					// Current choice is the clear winner
					continue;
				}

				// Compute the optimal actions by comparing the intervals
				allowedChoices = getAllowedActions(state);
				while (allowedChoices.hasNext()) {
					int choice = allowedChoices.nextInt();
					if (comparator.compare(choiceBestBounds.get(choice), otherOptimalWorstBound) >= 0) {
						stateOptimalActions.set(choice);
					}
				}

				if (comparator.compare(currentChoiceBestBound, otherOptimalWorstBound) < 0) {
					// Current choice is a clear loser
					improvedStates.set(state);
					assert otherOptimalBestBoundChoice >= 0;
					policy.setChoiceIndex(state, otherOptimalBestBoundChoice);
					improvements++;
					continue;
				}
				// Results are not conclusive, the intervals are overlapping
				assert min && (Math.max(currentChoiceBestBound, otherOptimalBestBound) <= Math.min(currentChoiceWorstBound, otherOptimalBestBound))
						|| !min && (Math.max(currentChoiceWorstBound, otherOptimalWorstBound) <= Math.min(currentChoiceBestBound, otherOptimalWorstBound));
				unknownStates.set(state);
			}
			this.stateImprovements += improvements;
		}
	}
}
