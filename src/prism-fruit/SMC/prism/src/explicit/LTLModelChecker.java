//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Joachim Klein <klein@tcs.inf.tu-dresden.de> (TU Dresden)
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

import acceptance.AcceptanceBuchi;
import acceptance.AcceptanceGenRabin;
import acceptance.AcceptanceOmega;
import acceptance.AcceptanceRabin;
import acceptance.AcceptanceStreett;
import acceptance.AcceptanceType;
import automata.DA;
import automata.LTL2DA;
import common.FastUtils;
import common.Time;
import de.tum.in.naturals.set.BoundedNatBitSet;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import owl.LDBAWrapper;
import owl.automaton.edge.Edge;
import parser.State;
import parser.VarList;
import parser.ast.Declaration;
import parser.ast.DeclarationInt;
import parser.ast.Expression;
import parser.ast.ExpressionBinaryOp;
import parser.ast.ExpressionLabel;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.type.TypeBool;
import parser.type.TypePathBool;
import prism.ModelType;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismNotSupportedException;
import prism.PrismUtils;

import java.awt.Point;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.function.IntConsumer;

/**
 * LTL model checking functionality
 */
public class LTLModelChecker extends PrismComponent
{
	/** Make LTL product accessible as a Product */
	public class LTLProduct<M extends Model> extends Product<M>
	{
		private int daSize;
		private int invMap[];
		private AcceptanceOmega acceptance;

		public LTLProduct(M productModel, M originalModel, AcceptanceOmega acceptance, int daSize, int[] invMap)
		{
			super(productModel, originalModel);
			this.daSize = daSize;
			this.invMap = invMap;
			this.acceptance = acceptance;
		}

		@Override
		public int getModelState(int productState)
		{
			return invMap[productState] / daSize;
		}

		@Override
		public int getAutomatonState(int productState)
		{
			return invMap[productState] % daSize;
		}

		public AcceptanceOmega getAcceptance() {
			return acceptance;
		}

		public void setAcceptance(AcceptanceOmega acceptance) {
			this.acceptance = acceptance;
		}
	}

	/**
	 * Make LTL product accessible as a Product
	 *
	 * @author sickert
	 */
	public class LDBAProduct<M extends NondetModel> extends Product<M>
	{
		private final Int2IntMap invMap;
		private final NatBitSet goalStates;

		public LDBAProduct(M productModel, M originalModel, Int2IntMap invMap, NatBitSet goalStates)
		{
			super(productModel, originalModel);
			this.invMap = invMap;
			invMap.defaultReturnValue(-1);
			this.goalStates = goalStates;
		}

		@Override
		public int getAutomatonState(int productState)
		{
			return invMap.get(productState);
		}

		public NatBitSet getGoalStates()
		{
			return goalStates;
		}

		@Override
		public int getModelState(int productState)
		{
			return invMap.get(productState);
		}
	}

	/**
	 * @author sickert
	 */
	private final class LDBAProductBuilder {
		private final LDBAWrapper da;
		private final MDP model;
		private final int daSize;
		private final List<NatBitSet> labelBS;
		private final NatBitSet statesOfInterest;

		// Encoding:
		// each state s' = <s, q> = s + (modelSize * q)
		// s(s') = s' % modelSize
		// q(s') = s' / modelSize
		private int[] map;
		private int trapState;

		private LDBAProductBuilder(MDP m, LDBAWrapper da, List<NatBitSet> labelBS, NatBitSet statesOfInterest) {
			model = m;
			this.da = da;
			daSize = da.size();

			this.labelBS = labelBS;
			this.statesOfInterest = statesOfInterest;
			this.map = new int[m.getNumStates() * da.size()];
			Arrays.fill(map, -1);
			trapState = -1;
		}

		/**
		 * Lookup the product state
		 * @param s The model state
		 * @param q The LDBA state
		 */
		private int getIndex(int s, int q) {
			return (s * daSize) + q;
		}

		private void buildTrapState(MDPSimple prodModel) {
			trapState = prodModel.addState();
			prodModel.addChoice(trapState, new Distribution(trapState, 1.0d));
		}

		private LDBAProduct<MDP> build() throws PrismException {
			MDPSimple prodModel = new MDPSimple();

			if (model.getVarList() != null) {
				VarList varList = model.getVarList();
				// Create a (new, unique) name for the variable that will represent DA states
				String daVar = "_da";
				while (varList.getIndex(daVar) != -1) {
					daVar = "_" + daVar;
				}

				VarList newVarList = (VarList) varList.clone();
				Declaration decl = new Declaration(daVar, new DeclarationInt(Expression.Int(0), Expression.Int(Math.max(da.size() - 1, 1))));
				newVarList.addVar(0, decl, 1, model.getConstantValues());
				prodModel.setVarList(newVarList);
			}

			Deque<Point> queue = new ArrayDeque<>();

			// We need results for all states of the original model in statesOfInterest
			// We thus explore states of the product starting from these states.
			// These are designated as initial states of the product model
			// (a) to ensure reachability is done for these states; and
			// (b) to later identify the corresponding product state for the original states
			//     of interest

			Int2IntMap initialStateInverseMapping = new Int2IntArrayMap(statesOfInterest.size());

			IntIterator iterator = statesOfInterest.iterator();
			while (iterator.hasNext()) {
				int s0 = iterator.nextInt();
				// Find corresponding initial state in DA
				int q0 = da.getEdgeDestByLabel(LDBAWrapper.getInitialState(), getLabel(s0));

				if (q0 >= 0) {
					int pi = getIndex(s0, q0);
					int p0 = map[pi];

					if (p0 == -1) {
						p0 = prodModel.addState();
						map[pi] = p0;
						queue.add(new Point(s0, q0));
						prodModel.addInitialState(p0);
						initialStateInverseMapping.put(p0, s0);
					}
					// TODO: detect "true" trap state. -> stop constructing and mark as goal state.
				} else if (trapState == -1) {
					buildTrapState(prodModel);
					prodModel.addInitialState(trapState);
					initialStateInverseMapping.put(trapState, s0);
				}
			}

			boolean initialComponentPresent = (da.getInitialComponentSize() > 0);
			NatBitSet goalStates = NatBitSets.set();
			constructReachableStates(prodModel, queue, goalStates);

			if (initialComponentPresent) {
				// Build a mapping from state indices to states (s,q), encoded as (s + (modelSize * q))
				int[] invMap = new int[prodModel.getNumStates()];
				for (int i1 = da.getInitialComponentSize() - 1; i1 >= 0; i1--) {
					int j = map[i1];

					if (j != -1) {
						invMap[j] = i1;
					}
				}

				// MECs have to be disjoint!
				for (MEC mec : ECComputerFast.computeMECs(prodModel)) {
					mec.states.forEach((IntConsumer) i -> {
						// Get Model state
						int s0 = invMap[i] / daSize;
						// Get LDBA state.
						int q0 = invMap[i] % daSize;

						int[] js = da.getJumpTargets(q0);

						for (int j : js) {
							int pi = getIndex(s0, j);
							int p = map[pi];

							// Add state to model
							if (p == -1) {
								p = prodModel.addState();
								map[pi] = p;
								queue.add(new Point(s0, j));
								constructReachableStates(prodModel, queue, goalStates);
							}

							// Add transition
							prodModel.addChoice(i, new Distribution(p, 1.0d));
						}
					});
				}
			}

			// Clear memory
			map = null;

			// Copy model to the sparse model to improve performance
			return new LDBAProduct<>(new MDPSparse(prodModel), model, initialStateInverseMapping, goalStates);
		}

		// Get NatBitSet representing APs (labels) satisfied by successor state s
		private NatBitSet getLabel(int s) {
			NatBitSet label = NatBitSets.boundedSet(labelBS.size());
			ListIterator<NatBitSet> iterator = labelBS.listIterator();

			while (iterator.hasNext()) {
				if (iterator.next().contains(s)) {
					label.set(iterator.previousIndex());
				}
			}

			return label;
		}

		private void constructReachableStates(MDPSimple prodModel, Deque<Point> queue, NatBitSet goalStates) {
			// local acceptance
			List<Int2ObjectMap<NatBitSet>> buchi = new ArrayList<>();
			NatBitSet preMEC = null;

			if (goalStates != null) {
				preMEC = NatBitSets.set();
			}

			for (int i = 0; i < da.getAcceptanceSize(); i++) {
				buchi.add(new Int2ObjectLinkedOpenHashMap<>());
			}

			while (!queue.isEmpty()) {
				Point current = queue.remove();

				int s1 = current.x;
				int q1 = current.y;
				int p1 = map[getIndex(s1, q1)];

				if (preMEC != null) {
					preMEC.set(p1);
				}

				// Go through transitions from state s_1 in original model
				final int numChoices = model.getNumChoices(s1);

				for (int j = 0; j < numChoices; j++) {
					Distribution prodDistr = new Distribution();
					NatBitSet accIndexes = NatBitSets.set();

					for (Int2DoubleMap.Entry e : model.getTransitions(s1, j)) {
						int s2 = e.getIntKey();
						NatBitSet label = getLabel(s2);
						Edge<?> edge = da.getEdgeByLabel(q1, label);

						int p2;

						if (edge != null) {
							int q2 = da.get(edge.getSuccessor());
							int pi = getIndex(s2, q2);
							p2 = map[pi];

							// Add state to model
							if (p2 == -1) {
								p2 = prodModel.addState();
								map[pi] = p2;
								queue.add(new Point(s2, q2));
							}

							// If there is an acceptance mapping, select all accepting edges
							edge.acceptanceSetIterator().forEachRemaining((IntConsumer) accIndexes::set);
						} else if (trapState == -1) {
							buildTrapState(prodModel);
							p2 = trapState;
						} else {
							p2 = trapState;
						}

						// Add transition
						prodDistr.add(p2, e.getDoubleValue());
					}

					IntIterator iterator = accIndexes.iterator();
					while (iterator.hasNext()) {
						int i = iterator.nextInt();
						Int2ObjectMap<NatBitSet> map = buchi.get(i);
						NatBitSet actions = map.computeIfAbsent(p1, k -> NatBitSets.boundedSet(numChoices));
						actions.set(j);
					}

					prodModel.addChoice(p1, prodDistr);
				}
			}

			if (goalStates != null && !preMEC.isEmpty()) {
				NatBitSet localGoalStates = findAcceptingECStatesForLDBA(prodModel, buchi, Collections.singleton(MEC.createMEC(prodModel, preMEC)));
				goalStates.or(localGoalStates);
			}
		}
	}

	/**
	 * Create a new LTLModelChecker, inherit basic state from parent (unless null).
	 */
	public LTLModelChecker(PrismComponent parent)
	{
		super(parent);
	}

	/**
	 * Returns {@code true} if expression {@code expr} is a formula that can be handled by
	 * LTLModelChecker for the given ModelType.
	 */
	public static boolean isSupportedLTLFormula(ModelType modelType, Expression expr) throws PrismLangException
	{
		if (!expr.isPathFormula(true)) {
			return false;
		}
		if (Expression.containsTemporalTimeBounds(expr)) {
			if (modelType.continuousTime()) {
				// Only support temporal bounds for discrete time models
				return false;
			}

			if (!expr.isSimplePathFormula()) {
				// Only support temporal bounds for simple path formulas
				return false;
			}
		}
		return true;
	}

	/**
	 * Extract maximal state formula from an LTL path formula, model check them (with passed in model checker) and
	 * replace them with ExpressionLabel objects L0, L1, etc. Expression passed in is modified directly, but the result
	 * is also returned. As an optimisation, expressions that results in true/false for all states are converted to an
	 * actual true/false, and duplicate results (or their negations) reuse the same label. NatBitSets giving the states which
	 * satisfy each label are put into the vector labelBS, which should be empty when this function is called.
	 */
	public Expression checkMaximalStateFormulas(StateModelChecker mc, Model model, Expression expr, List<NatBitSet> labelBS) throws PrismException
	{
		// A state formula
		if (expr.getType() instanceof TypeBool) {
			// Model check state formula for all states
			StateValues sv = mc.checkExpression(model, expr, null);
			NatBitSet bs = sv.getNatBitSet();
			// Detect special cases (true, false) for optimisation
			if (bs.isEmpty()) {
				return Expression.False();
			}
			if (bs.size() == model.getNumStates()) {
				return Expression.True();
			}
			// See if we already have an identical result
			// (in which case, reuse it)
			int i = labelBS.indexOf(bs);
			if (i != -1) {
				sv.clear();
				return new ExpressionLabel("L" + i);
			}
			// Also, see if we already have the negation of this result
			// (in which case, reuse it)
			NatBitSet bsNeg = NatBitSets.boundedFilledSet(model.getNumStates());
			bsNeg.andNot(bs);
			i = labelBS.indexOf(bsNeg);
			if (i != -1) {
				sv.clear();
				return Expression.Not(new ExpressionLabel("L" + i));
			}
			// Otherwise, add result to list, return new label
			labelBS.add(bs);
			return new ExpressionLabel("L" + (labelBS.size() - 1));
		}
		// A path formula (recurse, modify, return)
		if (expr.getType() instanceof TypePathBool) {
			if (expr instanceof ExpressionBinaryOp) {
				ExpressionBinaryOp exprBinOp = (ExpressionBinaryOp) expr;
				exprBinOp.setOperand1(checkMaximalStateFormulas(mc, model, exprBinOp.getOperand1(), labelBS));
				exprBinOp.setOperand2(checkMaximalStateFormulas(mc, model, exprBinOp.getOperand2(), labelBS));
			} else if (expr instanceof ExpressionUnaryOp) {
				ExpressionUnaryOp exprUnOp = (ExpressionUnaryOp) expr;
				exprUnOp.setOperand(checkMaximalStateFormulas(mc, model, exprUnOp.getOperand(), labelBS));
			} else if (expr instanceof ExpressionTemporal) {
				ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
				if (exprTemp.getOperand1() != null) {
					exprTemp.setOperand1(checkMaximalStateFormulas(mc, model, exprTemp.getOperand1(), labelBS));
				}
				if (exprTemp.getOperand2() != null) {
					exprTemp.setOperand2(checkMaximalStateFormulas(mc, model, exprTemp.getOperand2(), labelBS));
				}
			}
		}
		return expr;
	}

	/**
	 * Construct a deterministic automaton (DA) for an LTL formula, having first extracted maximal state formulas
	 * and model checked them with the passed in model checker. The maximal state formulas are assigned labels
	 * (L0, L1, etc.) which become the atomic propositions in the resulting DA. NatBitSets giving the states which
	 * satisfy each label are put into the vector {@code labelBS}, which should be empty when this function is called.
	 *
	 * @param mc a ProbModelChecker, used for checking maximal state formulas
	 * @param model the model
	 * @param expr a path expression, i.e. the LTL formula
	 * @param labelBS empty vector to be filled with NatBitSets for subformulas
	 * @param allowedAcceptance the allowed acceptance types
	 * @return the DA
	 */
	public DA<NatBitSet, ? extends AcceptanceOmega> constructDAForLTLFormula(ProbModelChecker mc, Model model, Expression expr, List<NatBitSet> labelBS,
			AcceptanceType... allowedAcceptance) throws PrismException
	{
		if (Expression.containsTemporalTimeBounds(expr)) {
			if (model.getModelType().continuousTime()) {
				throw new PrismException("Automaton construction for time-bounded operators not supported for " + model.getModelType()+".");
			}

			if (!expr.isSimplePathFormula()) {
				throw new PrismNotSupportedException("Time-bounded operators not supported in LTL: " + expr);
			}
		}

		// Model check maximal state formulas
		Expression ltl = checkMaximalStateFormulas(mc, model, expr.deepCopy(), labelBS);

		// Convert LTL formula to deterministic automaton
		mainLog.println("\nBuilding deterministic automaton (for " + ltl + ")...");
		long time = System.currentTimeMillis();
		LTL2DA ltl2da = new LTL2DA(this);
		DA<NatBitSet, ? extends AcceptanceOmega> da = ltl2da.convertLTLFormulaToDA(ltl, mc.getConstantValues(), allowedAcceptance);
		mainLog.println(da.getAutomataType() + " has " + da.size() + " states, " + da.getAcceptance().getSizeStatistics() + ".");
		da.checkForCanonicalAPs(labelBS.size());
		time = System.currentTimeMillis() - time;
		mainLog.println("Time for "+da.getAutomataType()+" translation: " + time / 1000.0 + " seconds.");
		// If required, export DA
		if (settings.getExportPropAut()) {
			mainLog.println("Exporting " + da.getAutomataType() + " to file \"" + settings.getExportPropAutFilename() + "\"...");
			PrintStream out = PrismUtils.newPrintStream(settings.getExportPropAutFilename());
			da.print(out, settings.getExportPropAutType());
			out.close();
		}

		return da;
	}

	/**
	 * Construct an LDBA for an LTL formula, having first extracted maximal state formulas
	 * and model checked them with the passed in model checker. The maximal state formulas are assigned labels
	 * (L0, L1, etc.) which become the atomic propositions in the resulting DA. NatBitSets giving the states which
	 * satisfy each label are put into the vector {@code labelBS}, which should be empty when this function is called.
	 *
	 * @param mc a ProbModelChecker, used for checking maximal state formulas
	 * @param model the model
	 * @param expr a path expression, i.e. the LTL formula
	 * @param labelBS empty vector to be filled with NatBitSets for subformulas
	 * @return the LDBA
	 * @author sickert
	 */
	private LDBAWrapper constructLDBAForLTLFormula(ProbModelChecker mc, MDP model, Expression expr, List<NatBitSet> labelBS) throws PrismException
	{
		if (Expression.containsTemporalTimeBounds(expr)) {
			if (model.getModelType().continuousTime()) {
				throw new PrismException("Automaton construction for time-bounded operators not supported for " + model.getModelType()+".");
			}

			if (!expr.isSimplePathFormula()) {
				throw new PrismNotSupportedException("Time-bounded operators not supported in LTL: " + expr);
			}
		}
		return LDBAWrapper.create(checkMaximalStateFormulas(mc, model, expr.deepCopy(), labelBS), mainLog);
	}

	/**
	 * Generate a deterministic automaton for the given LTL formula
	 * and construct the product of this automaton with a Markov chain.
	 *
	 * @param mc a ProbModelChecker, used for checking maximal state formulas
	 * @param model the model
	 * @param expr a path expression
 	 * @param statesOfInterest the set of states for which values should be calculated (null = all states)
 	 * @param allowedAcceptance the allowed acceptance types
	 * @return the product with the DA
	 */
	public LTLProduct<DTMC> constructProductMC(ProbModelChecker mc, DTMC model, Expression expr, NatBitSet statesOfInterest, AcceptanceType... allowedAcceptance) throws PrismException
	{
		// Convert LTL formula to automaton
		List<NatBitSet> labelBS = new ArrayList<>();
		DA<NatBitSet,? extends AcceptanceOmega> da;
		da = constructDAForLTLFormula(mc, model, expr, labelBS, allowedAcceptance);

		// Build product of model and automaton
		mainLog.println("\nConstructing MC-"+da.getAutomataType()+" product...");
		LTLProduct<DTMC> product = constructProductModel(da, model, labelBS, statesOfInterest);
		mainLog.print("\n" + product.getProductModel().infoStringTable());

		return product;
	}

	/**
	 * Generate a deterministic automaton for the given LTL formula
	 * and construct the product of this automaton with an MDP.
	 *
	 * @param mc a ProbModelChecker, used for checking maximal state formulas
	 * @param model the model
	 * @param expr a path expression
	 * @param statesOfInterest the set of states for which values should be calculated (null = all states)
	 * @param allowedAcceptance the allowed acceptance conditions
	 * @return the product with the DA
	 * @throws PrismException
	 */
	public LTLProduct<MDP> constructProductMDP(ProbModelChecker mc, MDP model, Expression expr, NatBitSet statesOfInterest, AcceptanceType... allowedAcceptance) throws PrismException
	{
		// Convert LTL formula to automaton
		List<NatBitSet> labelBS = new ArrayList<>();
		DA<NatBitSet,? extends AcceptanceOmega> da = constructDAForLTLFormula(mc, model, expr, labelBS, allowedAcceptance);

		// Build product of model and automaton
		Time time = new Time();
		mainLog.println("\nConstructing MDP-"+da.getAutomataType()+" product...");
		LTLProduct<MDP> product = constructProductModel(da, model, labelBS, statesOfInterest);
		double elapsedSeconds = time.elapsedSeconds();
		mainLog.print(String.format("Construction Time: %f seconds.%n", elapsedSeconds));
		mainLog.print("\n" + product.getProductModel().infoStringTable());

		return product;
	}

	/**
	 * Generate a limit-deterministic buchi automaton (LDBA) for the given LTL formula
	 * and construct the product of this automaton with an MDP.
	 *
	 * @param mc a ProbModelChecker, used for checking maximal state formulas
	 * @param model the model
	 * @param expr a path expression
	 * @param statesOfInterest the set of states for which values should be calculated (null = all states)
	 * @return the product with the LDBA
	 * @author sickert
	 */
	public LDBAProduct<MDP> constructProductMDPLDBA(ProbModelChecker mc, MDP model, Expression expr, NatBitSet statesOfInterest) throws PrismException
	{
		// Convert LTL formula to automaton
		List<NatBitSet> labelBS = new ArrayList<>();
		LDBAWrapper da = constructLDBAForLTLFormula(mc, model, expr, labelBS);

		// Build product of model and automaton
		Time time = new Time();
		mainLog.println("\nConstructing MDP-LDBA product and finding accepting MECs...");
		LDBAProductBuilder builder = new LDBAProductBuilder(model, da, labelBS, statesOfInterest);
		LDBAProduct<MDP> product = builder.build();
		double elapsedSeconds = time.elapsedSeconds();
		mainLog.print(String.format("Construction Time: %f seconds.%n", elapsedSeconds));
		mainLog.print("\n" + product.getProductModel().infoStringTable());

		return product;
	}

	/**
	 * Generate a deterministic automaton for the given LTL formula
	 * and construct the product of this automaton with an STPG.
	 *
	 * @param mc a ProbModelChecker, used for checking maximal state formulas
	 * @param model the model
	 * @param expr a path expression
	 * @param statesOfInterest the set of states for which values should be calculated (null = all states)
	 * @param allowedAcceptance the allowed acceptance conditions
	 * @return the product with the DA
	 * @throws PrismException
	 */
	public LTLProduct<STPG> constructProductSTPG(ProbModelChecker mc, STPG model, Expression expr, NatBitSet statesOfInterest, AcceptanceType... allowedAcceptance) throws PrismException
	{
		// Convert LTL formula to automaton
		List<NatBitSet> labelBS = new ArrayList<>();
		DA<NatBitSet,? extends AcceptanceOmega> da;
		da = constructDAForLTLFormula(mc, model, expr, labelBS, allowedAcceptance);

		// Build product of model and automaton
		mainLog.println("\nConstructing STPG-"+da.getAutomataType()+" product...");
		LTLProduct<STPG> product = constructProductModel(da, model, labelBS, statesOfInterest);
		mainLog.print("\n" + product.getProductModel().infoStringTable());

		return product;
	}

	/**
	 * Generate a deterministic automaton for the given LTL formula
	 * and construct the product of this automaton with a model.
	 *
	 * @param mc a ProbModelChecker, used for checking maximal state formulas
	 * @param model the model
	 * @param expr a path expression
	 * @param statesOfInterest the set of states for which values should be calculated (null = all states)
	 * @param allowedAcceptance the allowed acceptance conditions
	 * @return the product with the DA
	 * @throws PrismException
	 */
	public <M extends Model> LTLProduct<M> constructProductModel(ProbModelChecker mc, M model, Expression expr, NatBitSet statesOfInterest, AcceptanceType... allowedAcceptance) throws PrismException
	{
		// Convert LTL formula to automaton
		List<NatBitSet> labelBS = new ArrayList<>();
		DA<NatBitSet,? extends AcceptanceOmega> da;
		da = constructDAForLTLFormula(mc, model, expr, labelBS, allowedAcceptance);

		// Build product of model and automaton
		mainLog.println("\nConstructing " + model.getModelType() + "-" + da.getAutomataType() + " product...");
		LTLProduct<M> product = constructProductModel(da, model, labelBS, statesOfInterest);
		mainLog.print("\n" + product.getProductModel().infoStringTable());

		return product;
	}

	/**
	 * Construct the product of a DA and a model.
	 * @param da The DA
	 * @param model The model
	 * @param labelBS NatBitSets giving the set of states for each AP in the DA
	 * @param statesOfInterest the set of states for which values should be calculated (null = all states)
	 * @return The product model
	 */
	public <M extends Model> LTLProduct<M> constructProductModel(DA<NatBitSet, ? extends AcceptanceOmega> da, M model, List<NatBitSet> labelBS,
			NatBitSet statesOfInterest) throws PrismException
	{
		ModelType modelType = model.getModelType();
		int daSize = da.size();
		int numAPs = da.getAPList().size();
		int modelNumStates = model.getNumStates();
		int prodNumStates = modelNumStates * daSize;
		int s_1, s_2, q_1, q_2;
		NatBitSet s_labels = NatBitSets.boundedSet(numAPs);
		List<State> prodStatesList = null, daStatesList = null;

		VarList newVarList = null;

		if (model.getVarList() != null) {
			VarList varList = model.getVarList();
			// Create a (new, unique) name for the variable that will represent DA states
			String daVar = "_da";
			while (varList.getIndex(daVar) != -1) {
				daVar = "_" + daVar;
			}

			newVarList = (VarList) varList.clone();
			// NB: if DA only has one state, we add an extra dummy state
			Declaration decl = new Declaration(daVar, new DeclarationInt(Expression.Int(0), Expression.Int(Math.max(da.size() - 1, 1))));
			newVarList.addVar(0, decl, 1, model.getConstantValues());
		}

		// Create a (simple, mutable) model of the appropriate type
		ModelSimple prodModel = null;
		switch (modelType) {
		case DTMC: {
			DTMCSimple dtmcProd = new DTMCSimple();
			dtmcProd.setVarList(newVarList);
			prodModel = dtmcProd;
			break;
		}
		case MDP: {
			MDPSimple mdpProd = new MDPSimple();
			mdpProd.setVarList(newVarList);
			prodModel = mdpProd;
			break;
		}
		case STPG: {
			STPGExplicit stpgProd = new STPGExplicit();
			stpgProd.setVarList(newVarList);
			prodModel = stpgProd;
			break;
		}
		default:
			throw new PrismNotSupportedException("Model construction not supported for " + modelType + "s");
		}

		// Encoding:
		// each state s' = <s, q> = s * daSize + q
		// s(s') = s' / daSize
		// q(s') = s' % daSize

		LinkedList<Point> queue = new LinkedList<>();
		int map[] = new int[prodNumStates];
		Arrays.fill(map, -1);

		if (model.getStatesList() != null) {
			prodStatesList = new ArrayList<>();
			daStatesList = new ArrayList<>(da.size());
			for (int i = 0; i < da.size(); i++) {
				daStatesList.add(new State(1).setValue(0, i));
			}
		}

		// We need results for all states of the original model in statesOfInterest
		// We thus explore states of the product starting from these states.
		// These are designated as initial states of the product model
		// (a) to ensure reachability is done for these states; and
		// (b) to later identify the corresponding product state for the original states
		//     of interest
		IntIterator iterator = FastUtils.iterator(statesOfInterest, model.getNumStates());
		while (iterator.hasNext()) {
			int s_0 = iterator.nextInt();
			// Get NatBitSet representing APs (labels) satisfied by state s_0
			for (int k = 0; k < numAPs; k++) {
				s_labels.set(k, labelBS.get(Integer.parseInt(da.getAPList().get(k).substring(1))).contains(s_0));
			}
			// Find corresponding initial state in DA
			int q_0 = da.getEdgeDestByLabel(da.getStartState(), s_labels);
			if (q_0 < 0) {
				throw new PrismException("The deterministic automaton is not complete (state " + da.getStartState() + ")");
			}
			// Add (initial) state to product
			queue.add(new Point(s_0, q_0));
			switch (modelType) {
			case STPG:
				((STPGExplicit) prodModel).addState(((STPG) model).getPlayer(s_0));
				break;
			default:
				prodModel.addState();
			break;
			}
			prodModel.addInitialState(prodModel.getNumStates() - 1);
			map[s_0 * daSize + q_0] = prodModel.getNumStates() - 1;
			if (prodStatesList != null) {
				// Store state information for the product
				prodStatesList.add(new State(daStatesList.get(q_0), model.getStatesList().get(s_0)));
			}
		}

		// Product states
		NatBitSet visited = NatBitSets.boundedSet(prodNumStates);
		while (!queue.isEmpty()) {
			Point p = queue.pop();
			s_1 = p.x;
			q_1 = p.y;
			visited.set(s_1 * daSize + q_1);

			// Go through transitions from state s_1 in original model
			int numChoices = (model instanceof NondetModel) ? ((NondetModel) model).getNumChoices(s_1) : 1;
			for (int j = 0; j < numChoices; j++) {
				Iterator<Int2DoubleMap.Entry> iter;
				switch (modelType) {
				case DTMC:
					iter = ((DTMC) model).getTransitionsIterator(s_1);
					break;
				case MDP:
					iter = ((MDP) model).getTransitions(s_1, j).iterator();
					break;
				case STPG:
					iter = ((STPG) model).getTransitions(s_1, j).iterator();
					break;
				default:
					throw new PrismNotSupportedException("Product construction not implemented for " + modelType + "s");
				}
				Distribution prodDistr = null;
				if (modelType.nondeterministic()) {
					prodDistr = new Distribution();
				}
				while (iter.hasNext()) {
					Int2DoubleMap.Entry e = iter.next();
					s_2 = e.getIntKey();
					double prob = e.getDoubleValue();
					// Get NatBitSet representing APs (labels) satisfied by successor state s_2
					for (int k = 0; k < numAPs; k++) {
						s_labels.set(k, labelBS.get(Integer.parseInt(da.getAPList().get(k).substring(1))).contains(s_2));
					}
					// Find corresponding successor in DA
					q_2 = da.getEdgeDestByLabel(q_1, s_labels);
					if (q_2 < 0) {
						throw new PrismException("The deterministic automaton is not complete (state " + q_1 + ")");
					}
					// Add state/transition to model
					if (!visited.contains(s_2 * daSize + q_2) && map[s_2 * daSize + q_2] == -1) {
						queue.add(new Point(s_2, q_2));
						switch (modelType) {
						case STPG:
							((STPGExplicit) prodModel).addState(((STPG) model).getPlayer(s_2));
							break;
						default:
							prodModel.addState();
							break;
						}
						map[s_2 * daSize + q_2] = prodModel.getNumStates() - 1;
						if (prodStatesList != null) {
							// Store state information for the product
							prodStatesList.add(new State(daStatesList.get(q_2), model.getStatesList().get(s_2)));
						}
					}
					switch (modelType) {
					case DTMC:
						((DTMCSimple) prodModel).setProbability(map[s_1 * daSize + q_1], map[s_2 * daSize + q_2], prob);
						break;
					case MDP:
					case STPG:
						prodDistr.set(map[s_2 * daSize + q_2], prob);
						break;
					default:
						throw new PrismNotSupportedException("Product construction not implemented for " + modelType + "s");
					}
				}
				switch (modelType) {
				case MDP:
					((MDPSimple) prodModel).addActionLabelledChoice(map[s_1 * daSize + q_1], prodDistr, ((MDP) model).getAction(s_1, j));
					break;
				case STPG:
					((STPGExplicit) prodModel).addActionLabelledChoice(map[s_1 * daSize + q_1], prodDistr, ((STPG) model).getAction(s_1, j));
					break;
				default:
					break;
				}
			}
		}

		// Build a mapping from state indices to states (s,q), encoded as (s * daSize + q)
		int invMap[] = new int[prodModel.getNumStates()];
		for (int i = 0; i < map.length; i++) {
			if (map[i] != -1) {
				invMap[map[i]] = i;
			}
		}

		prodModel.findDeadlocks(false);

		if (prodStatesList != null) {
			prodModel.setStatesList(prodStatesList);
		}

		@SuppressWarnings("unchecked")
		LTLProduct<M> product = new LTLProduct<>((M) prodModel, model, null, daSize, invMap);

		// generate acceptance for the product model by lifting
		product.setAcceptance(liftAcceptance(product, da.getAcceptance()));

		// lift the labels
		for (String label : model.getLabels()) {
			NatBitSet liftedLabel = product.liftFromModel(model.getLabelStates(label));
			prodModel.addLabel(label, liftedLabel);
		}

		return product;
	}

	/**
	 * Find the set of states that belong to accepting BSCCs in a model wrt an acceptance condition.
	 * @param model The model
	 * @param acceptance The acceptance condition
	 */
	public NatBitSet findAcceptingBSCCs(Model model, AcceptanceOmega acceptance) throws PrismException
	{
		// Compute bottom strongly connected components (BSCCs)
		SCCComputer sccComputer = SCCComputer.createSCCComputer(this, model);
		sccComputer.computeBSCCs();
		List<NatBitSet> bsccs = sccComputer.getBSCCs();

		NatBitSet result = NatBitSets.boundedSet(model.getNumStates());

		for (NatBitSet bscc : bsccs) {
			if (acceptance.isBSCCAccepting(bscc)) {
				// this BSCC is accepting
				result.or(bscc);
			}
		}

		return result;
	}

	/**
	 * Compute the set of states in end components of the model that are accepting
	 * with regard to the acceptance condition.
	 * @param model the model
	 * @param acceptance the acceptance condition
	 * @return NatBitSet with the set of states that are accepting
	 */
	public NatBitSet findAcceptingECStates(NondetModel model, AcceptanceOmega acceptance) throws PrismException
	{
		if (acceptance instanceof AcceptanceBuchi) {
			return findAcceptingECStatesForBuchi(model, (AcceptanceBuchi) acceptance);
		}
		if (acceptance instanceof AcceptanceRabin) {
			return findAcceptingECStatesForRabin(model, (AcceptanceRabin) acceptance);
		}
		if (acceptance instanceof AcceptanceStreett) {
			return findAcceptingECStatesForStreett(model, (AcceptanceStreett) acceptance);
		}
		if (acceptance instanceof AcceptanceGenRabin) {
			return findAcceptingECStatesForGeneralizedRabin(model, (AcceptanceGenRabin) acceptance);
		}
		throw new PrismNotSupportedException("Computing end components for acceptance type '"+acceptance.getType()+"' currently not supported (explicit engine).");
	}

	/**
	 * Find the set of states in accepting end components (ECs) in a nondeterministic model wrt a Büchi acceptance condition.
	 * @param model The model
	 * @param acceptance The acceptance condition
	 */
	public NatBitSet findAcceptingECStatesForBuchi(NondetModel model, AcceptanceBuchi acceptance) throws PrismException
	{
		NatBitSet allAcceptingStates = NatBitSets.boundedSet(model.getNumStates());

		if (acceptance.getAcceptingStates().isEmpty()) {
			return allAcceptingStates;
		}

		// Compute accepting maximum end components (MECs)
		ECComputer ecComputer = ECComputer.createECComputer(this, model);
		ecComputer.computeMECStates();
		List<NatBitSet> mecs = ecComputer.getMECStates();
		// Union of accepting MEC states
		for (NatBitSet mec : mecs) {
			if (mec.intersects(acceptance.getAcceptingStates())) {
				allAcceptingStates.or(mec);
			}
		}

		return allAcceptingStates;
	}

	/**
	 * Find the set of states in accepting end components (ECs) in a nondeterministic model wrt a Rabin acceptance condition.
	 * @param model The model
	 * @param acceptance The acceptance condition
	 */
	public NatBitSet findAcceptingECStatesForRabin(NondetModel model, AcceptanceRabin acceptance) throws PrismException
	{
		NatBitSet allAcceptingStates = NatBitSets.boundedSet(model.getNumStates());
		int getNumStates = model.getNumStates();

		// Go through the DRA acceptance pairs (L_i, K_i)
		for (AcceptanceRabin.RabinPair pair : acceptance) {
			// Find model states *not* satisfying L_i
			NatBitSet bitsetLi = pair.getL();
			BoundedNatBitSet statesLi_not = NatBitSets.boundedSet(model.getNumStates());
			statesLi_not.orNot(bitsetLi);

			// Skip pairs with empty !L_i
			if (statesLi_not.isEmpty())
				continue;
			// Compute accepting maximum end components (MECs) in !L_i
			ECComputer ecComputer = ECComputer.createECComputer(this, model);
			ecComputer.computeMECStates(statesLi_not, pair.getK());
			List<NatBitSet> mecs = ecComputer.getMECStates();
			// Union MEC states
			for (NatBitSet mec : mecs) {
				allAcceptingStates.or(mec);
			}
		}

		return allAcceptingStates;
	}

	/**
	 * Find the set of states in accepting end components (ECs) in a nondeterministic model wrt a Streett acceptance condition.
	 * @param model The model
	 * @param acceptance The Streett acceptance condition
	 */
	public NatBitSet findAcceptingECStatesForStreett(NondetModel model, AcceptanceStreett acceptance) throws PrismException
	{
		class ECandPairs {
			NatBitSet MEC;
			NatBitSet activePairs;
		}

		NatBitSet allAcceptingStates = NatBitSets.boundedSet(model.getNumStates());
		NatBitSet allPairs = NatBitSets.boundedFilledSet(acceptance.size());

		Deque<ECandPairs> todo = new ArrayDeque<>();
		ECComputer ecComputer = ECComputer.createECComputer(this, model);
		ecComputer.computeMECStates();
		for (NatBitSet mecs : ecComputer.getMECStates()) {
			ECandPairs ecp = new ECandPairs();
			ecp.MEC = mecs;
			ecp.activePairs = allPairs;
			todo.push(ecp);
		}

		while (!todo.isEmpty()) {
			ECandPairs ecp = todo.pop();
			NatBitSet newActivePairs = ecp.activePairs.clone();
			NatBitSet restrict = null;

			// check for acceptance
			boolean allAccepting = true;
			IntIterator pairIterator = ecp.activePairs.iterator();
			while (pairIterator.hasNext()) {
				int pair = pairIterator.nextInt();

				if (!acceptance.get(pair).isBSCCAccepting(ecp.MEC)) {
					// this pair is not accepting
					if (restrict == null) {
						restrict = ecp.MEC.clone();
					}
					restrict.andNot(acceptance.get(pair).getR());
					newActivePairs.clear(pair);
					allAccepting = false;
				}
			}

			if (allAccepting) {
				allAcceptingStates.or(ecp.MEC);
			} else if (restrict.isEmpty()) {
				// nothing to do
			} else {
				ecComputer = ECComputer.createECComputer(this, model);
				ecComputer.computeMECStates(restrict);
				for (NatBitSet mecs : ecComputer.getMECStates()) {
					ECandPairs newEcp = new ECandPairs();
					newEcp.MEC = mecs;
					newEcp.activePairs = newActivePairs;
					todo.push(newEcp);
				}
			}
		}

		return allAcceptingStates;
	}

	/**
	 * Find the set of states in accepting end components (ECs) in a nondeterministic model wrt a Generalized Rabin acceptance condition.
	 * @param model The model
	 * @param acceptance The acceptance condition
	 */
	public NatBitSet findAcceptingECStatesForGeneralizedRabin(NondetModel model, AcceptanceGenRabin acceptance) throws PrismException
	{
		NatBitSet allAcceptingStates = NatBitSets.boundedSet(model.getNumStates());
		int getNumStates = model.getNumStates();

		// Go through the GR acceptance pairs (L_i, K_i_1, ..., K_i_n)
		for (AcceptanceGenRabin.GenRabinPair pair : acceptance) {
			// Find model states *not* satisfying L_i
			NatBitSet bitsetLi = pair.getL();
			BoundedNatBitSet statesLi_not = NatBitSets.boundedSet(model.getNumStates());
			statesLi_not.orNot(bitsetLi);
			// Skip pairs with empty !L_i
			if (statesLi_not.isEmpty())
				continue;
			// Compute maximum end components (MECs) in !L_i
			ECComputer ecComputer = ECComputer.createECComputer(this, model);
			ecComputer.computeMECStates(statesLi_not);
			List<NatBitSet> mecs = ecComputer.getMECStates();
			// Check which MECs contain a state from each K_i_j
			int n = pair.getNumK();
			for (NatBitSet mec : mecs) {
				boolean allj = true;
				for (int j = 0; j < n; j++) {
					if (!mec.intersects(pair.getK(j))) {
						allj = false;
						break;
					}
				}
				if (allj) {
					allAcceptingStates.or(mec);
				}
			}
		}

		return allAcceptingStates;
	}

	/**
	 * Use the given Büchi condition to identify goal states in the model. Restrict the search to the given preMECs.
	 *
	 * @param model      The model MDP
	 * @param acceptance Generalised Büchi Acceptance
	 * @param preMECs    MEC candidates that need to be refined.
	 * @return a set of states needed to be reached for acceptance
	 * @author sickert
	 */
	public static NatBitSet findAcceptingECStatesForLDBA(NondetModel model, List<Int2ObjectMap<NatBitSet>> acceptance, Collection<MEC> preMECs)
	{
		NatBitSet allAcceptingStates = NatBitSets.boundedSet(model.getNumStates());

		for (MEC mec : ECComputerFast.computeMECs(model, preMECs)) {
			boolean allSuccess = true;

			for (Int2ObjectMap<NatBitSet> map : acceptance) {
				boolean indexSuccess = false;

				IntIterator iterator = mec.states.iterator();
				while (iterator.hasNext() && !indexSuccess) {
					int i = iterator.nextInt();
					NatBitSet actions = map.get(i);
					if (actions != null) {
						indexSuccess = actions.intersects(mec.actions.get(i));
					}
				}

				if (!indexSuccess) {
					allSuccess = false;
					break;
				}
			}

			if (allSuccess) {
				allAcceptingStates.or(mec.states);
			}
		}

		return allAcceptingStates;
	}

	/** Lift the acceptance condition from the automaton to the product states. */
	private AcceptanceOmega liftAcceptance(LTLProduct<?> product, AcceptanceOmega acceptance)
	{
		// make a copy of the acceptance condition
		AcceptanceOmega lifted = acceptance.clone();

		// lift state sets
		lifted.lift(new AcceptanceOmega.LiftNatBitSet() {
			@Override
			public NatBitSet lift(NatBitSet states)
			{
				return product.liftFromAutomaton(states);
			}
		});

		return lifted;
	}

}
