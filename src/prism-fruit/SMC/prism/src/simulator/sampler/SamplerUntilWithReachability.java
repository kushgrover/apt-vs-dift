package simulator.sampler;

import cern.colt.Timer;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.StateValues;
import it.unimi.dsi.fastutil.ints.IntArrayFIFOQueue;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntPriorityQueue;
import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;
import parser.State;
import parser.VarList;
import parser.ast.Expression;
import parser.ast.ExpressionExists;
import parser.ast.ExpressionLiteral;
import parser.ast.ExpressionTemporal;
import parser.ast.ModulesFile;
import prism.Model;
import prism.NonProbModelChecker;
import prism.Prism;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismLog;
import prism.PrismPrintStreamLog;
import simulator.Path;
import simulator.TransitionList;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.IntConsumer;

/**
 * Sampler for unbouded until properties that uses a reachability analysis to prune states that cannot satisfy the property.
 * Based on "Statistical Verification of Probabilistic Properties with Unbounded Until" Younes'11.
 *
 * @author Przemysław Daca
 */
public class SamplerUntilWithReachability extends SamplerBoolean
{
	// how often to check for states that cannot reach the target
	private final int bound = 100;
	private Expression left;
	private Expression right;
	private JDDNode reach;
	private Set<State> reachExpl;
	private boolean isExplicit;
	private int numVars;
	private VarList varList;
	// TODO check if JDD used correctly
	private JDDVars[] varDDRowVars;
	private Timer reachTmr;
	private Timer inclTmr;
	private int boundIdx;

	public SamplerUntilWithReachability(ModulesFile mf, ExpressionTemporal expr, int engine) throws PrismException
	{
		// Make sure expression is of the correct type
		// Then extract other required info
		if (expr.getOperator() != ExpressionTemporal.P_U) {
			throw new PrismException("Error creating Sampler");
		}
		left = expr.getOperand1();
		right = expr.getOperand2();

		reachTmr = new Timer();
		inclTmr = new Timer();

		reachTmr.start();
		// find states that satisfy [ E (true U right)]

		PrismLog mainLog = new PrismPrintStreamLog(System.out);
		Prism prism = new Prism(mainLog);
		prism.setEngine(engine);

		prism.loadPRISMModel(mf);
		prism.buildModel();
		ExpressionTemporal expFr = new ExpressionTemporal(ExpressionTemporal.P_U, ExpressionLiteral.True(), right);
		ExpressionExists expEFr = new ExpressionExists(expFr);

		if (engine == Prism.EXPLICIT) {
			explicit.Model model = prism.getBuiltModelExplicit();

			explicit.NonProbModelChecker mc = new explicit.NonProbModelChecker(prism);
			StateValues trgt = mc.checkExpression(model, right, null);
			NatBitSet trgtBs = trgt.getNatBitSet();

			// construct reverse tr. rel
			List<IntList> revtr = new ArrayList<>(model.getNumStates());
			for (int s = 0; s < model.getNumStates(); s++) {
				revtr.add(new IntArrayList());
			}

			for (int s = 0; s < model.getNumStates(); s++) {
				IntIterator iter = model.getSuccessorsIterator(s);
				while (iter.hasNext()) {
					int t = iter.nextInt();
					revtr.get(t).add(s);
				}
			}

			// do backward reachability from states that satisfy "right"
			NatBitSet done = NatBitSets.set();
			IntPriorityQueue worklist = new IntArrayFIFOQueue(trgtBs.size());
			trgtBs.forEach((IntConsumer) worklist::enqueue);

			while (!worklist.isEmpty()) {
				int s = worklist.dequeueInt();
				if (done.contains(s)) {
					continue;
				}
				done.set(s);
				revtr.get(s).forEach((IntConsumer) worklist::enqueue);
			}

			List<State> states = model.getStatesList();
			this.reachExpl = new HashSet<>();
			done.forEach((IntConsumer) s -> reachExpl.add(states.get(s)));
		} else {
			Model model = prism.getBuiltModel();
			numVars = model.getNumVars();
			varList = model.getVarList();
			varDDRowVars = model.getVarDDRowVars();

			NonProbModelChecker mc = new NonProbModelChecker(prism, model, null);

			// TODO remove this spurious command
			mc.check(ExpressionLiteral.True());
			reach = mc.checkExpressionDD(expEFr);
			//JDD.Ref(reach);
			model.getStartStates();
		}

		reachTmr.stop();

		// Initialise sampler info
		reset();
		resetStats();
	}

	@Override
	public void reset()
	{
		super.reset();
		boundIdx = 0;
	}

	@Override
	public void resetStats()
	{
		super.resetStats();
	}

	@Override
	public boolean update(Path path, TransitionList transList) throws PrismLangException
	{
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
		// Or, if we are now at a deadlock/self-loop
		else if (transList != null && (transList.isDeadlock() || path.isLooping())) {
			//else if (transList != null && (transList.isDeadlock() || transList.isDeterministicSelfLoop(currentState))) {
			valueKnown = true;
			value = false;
		} else {
			if (boundIdx++ == bound) {
				boundIdx = 0;
				if (!canReachRHS(currentState)) {
					// Stop if the state cannot reach RHS
					valueKnown = true;
					value = false;
				}
			}

		}
		// Otherwise, don't know

		return valueKnown;
	}

	/**
	 * True iff the RHS can be reached from the state.
	 *
	 * @param state
	 * @return
	 */
	private boolean canReachRHS(State state)
	{
		inclTmr.start();
		boolean res;
		if (isExplicit) {
			res = reachExpl.contains(state);
		} else {
			res = stateInDD(state, reach);
		}
		inclTmr.stop();
		return res;
	}

	/**
	 * Checks if the state belongs to the BDD.
	 *
	 * @param state
	 * @param dd
	 * @return
	 */
	private boolean stateInDD(State state, JDDNode dd)
	{
		int j = 0;
		//JDDNode res = JDD.Constant(1);

		JDDNode res = JDD.Constant(1);
		for (int i = 0; i < numVars; i++) {
			try {
				j = varList.encodeToInt(i, state.varValues[i]);
			} catch (PrismLangException e) {
				// Won't happen
			}
			res = JDD.Apply(JDD.TIMES, res, JDD.SetVectorElement(JDD.Constant(0), varDDRowVars[i], j, 1.0));
		}

		JDD.Ref(reach);
		res = JDD.And(reach, res);

		//JDD.ExportDDToDotFileLabelled(res, "state.dot", varNames)
		boolean b = !res.equals(JDD.ZERO);
		JDD.Deref(res);

		return b;
	}

	public String toString()
	{
		StringBuilder bldr = new StringBuilder();
		bldr.append(SamplerUntilWithReachability.class.getName() + " statistics:\n");
		bldr.append("Symbolic:                           " + !isExplicit + "\n");
		bldr.append("Time for reachability analysis:     " + reachTmr + "\n");
		bldr.append("Time for inclusion check:           " + inclTmr + "\n");
		return bldr.toString();
	}
}
