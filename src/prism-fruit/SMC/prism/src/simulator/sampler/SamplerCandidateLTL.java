package simulator.sampler;

import acceptance.AcceptanceOmega;
import automata.DA;
import automata.LTL2DA;
import cern.colt.Timer;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.LTLModelChecker;
import parser.State;
import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionBinaryOp;
import parser.ast.ExpressionLabel;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.ast.ModulesFile;
import parser.type.TypeBool;
import parser.type.TypePathBool;
import prism.Pair;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismUtils;
import simulator.Path;
import simulator.PathFull;
import simulator.TransitionList;
import simulator.UnboundedSimulationParameters;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Sampler checks if the trace satisfies an LTL property by identifying a BSCC candidate in the trace.
 *
 * @author Przemysław Daca
 */
public class SamplerCandidateLTL extends Sampler
{
	private static final boolean checkParam = true;
	// the initial visit number
	private final double baseVN;
	public Timer updateTmr = new Timer();
	public Timer cndTmr = new Timer();
	public Timer processTmr = new Timer();
	public Timer constructTmr = new Timer();
	public Timer checkTmr = new Timer();
	// Value of current path
	protected boolean value;
	// number of samples
	protected int numSamples;
	private Expression[] labelExp;
	// number of samples that satisfy the property
	private int numTrue;
	// the minimal visit number (VN) of status in the candidate
	private int visitNo;
	// how many times we increase VN after a new non-trivial candidate is found
	private double increseVN;
	// DRA for the property to be checked
	private DA<NatBitSet, ? extends AcceptanceOmega> dra;
	// last DRA state
	private int draState;
	// how often to check for a VN-candidate
	private int checkBound;
	private int checkIndx;
	// CHECKING CANDIDATES
	// last inspected index	(-1 if none)
	private int lastIndx;
	//private Candidate candidate;
	private Candidate<Pair<State, Integer>> lastNonTrivialCnd;
	// sequence of non-subsummed candidates  on the path
	private Deque<Candidate<Pair<State, Integer>>> cndLst;
	// state -> index in cndLst
	private Map<Pair<State, Integer>, Integer> cndMap;
	// STATISTICS
	// count the number of candidates that were false positives
	private boolean countFP = false;
	private long candLen;
	private long transLen;
	private long canTotal;
	private long transTotal;
	private long lastRunSteps;
	private int totalFP;
	private int entryToBirth;
	private int cndFound;
	private int cndFoundMax;
	private int cndFoundTotal;
	private int entryToBirthTotal;
	private double simMinProb;

	/**
	 * We need log base 2, not log base e.
	 * @param num Some rational number
	 * @return Ld(num), i.e. logarithm dualis
	 */
	public static double ld(double num){
		return Math.log(num)/Math.log(2);
	}

	public SamplerCandidateLTL(ModulesFile mf, Expression expr, Values constantValues, UnboundedSimulationParameters usp) throws PrismException
	{

		double simFalseCnd = usp.getSimFalseCnd();
		double ratio = usp.getSimRatio();
		this.checkBound = usp.getSimCheckBound();
		this.cndMap = new HashMap<>();
		this.cndLst = new ArrayDeque<>();

		//Maxi: I removed simFalseCnd from consideration, as it seemed like some weird optimization and I didn't understand it
		//Now I'm using it to decide between cautious and bold monitor
		this.simMinProb = usp.getSimMinProb();
		if (simFalseCnd > 0.5){
			//bold monitor, use sensible calculation for strongness of candidates (Recall visitNo basically is k_i)
			this.baseVN = ld(ratio) / ld(1 - simMinProb);
			this.increseVN = - (1 / Math.log(1 - simMinProb));
			visitNo = (int) Math.ceil(baseVN);
		}
		else{
			//cautious monitor: when having an x candidate, where x is fixed, check the value, and stop simulation if false.
			// if true, switch to checking whether to stop the experiment, using the bold formula (see update method)
			this.baseVN = 10*ratio;
			this.increseVN = 0;
			visitNo = (int) Math.ceil(baseVN);
		}

		if (Expression.containsTemporalTimeBounds(expr)) {
			throw new PrismException("Time-bounded operators not supported in LTL: " + expr);
		}

		List<Expression> labelBS = new ArrayList<>();
		LTLModelChecker ltlmc = new LTLModelChecker(null);
		ExpressionProb exprp = (ExpressionProb) expr;
		expr = checkMaximalStateFormulas(mf, exprp.getExpression().deepCopy(), labelBS);
		PrismComponent pc = new PrismComponent();
		LTL2DA ltl2da = new LTL2DA(pc);
		dra = ltl2da.convertLTLFormulaToDRA(expr, constantValues);
		this.labelExp = new Expression[labelBS.size()];
		labelBS.toArray(this.labelExp);
	}

	private Expression checkMaximalStateFormulas(ModulesFile mf, Expression expr, List<Expression> labelBS)
	{
		// A state formula
		if (expr.getType() instanceof TypeBool) {
			// TODO something smart, like detecting special cases
			labelBS.add(expr);
			return new ExpressionLabel("L" + (labelBS.size() - 1));
		}
		// A path formula (recurse, modify, return)
		else if (expr.getType() instanceof TypePathBool) {
			if (expr instanceof ExpressionBinaryOp) {
				ExpressionBinaryOp exprBinOp = (ExpressionBinaryOp) expr;
				exprBinOp.setOperand1(checkMaximalStateFormulas(mf, exprBinOp.getOperand1(), labelBS));
				exprBinOp.setOperand2(checkMaximalStateFormulas(mf, exprBinOp.getOperand2(), labelBS));
			} else if (expr instanceof ExpressionUnaryOp) {
				ExpressionUnaryOp exprUnOp = (ExpressionUnaryOp) expr;
				exprUnOp.setOperand(checkMaximalStateFormulas(mf, exprUnOp.getOperand(), labelBS));
			} else if (expr instanceof ExpressionTemporal) {
				ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
				if (exprTemp.getOperand1() != null) {
					exprTemp.setOperand1(checkMaximalStateFormulas(mf, exprTemp.getOperand1(), labelBS));
				}
				if (exprTemp.getOperand2() != null) {
					exprTemp.setOperand2(checkMaximalStateFormulas(mf, exprTemp.getOperand2(), labelBS));
				}
			}
		}
		return expr;
	}

	@Override
	public void reset()
	{
		valueKnown = false;
		value = false;
		//visitNo = (int) Math.ceil(baseVN);

		checkIndx = 0;
		lastIndx = -1;
		cndMap.clear();
		cndLst.clear();
		draState = dra.getStartState();
	}

	@Override
	public void resetStats()
	{
		numSamples = 0;
		numTrue = 0;
		candLen = 0;
		transLen = 0;
		//cndFound = 0;
		lastRunSteps = 0;
	}

	@Override
	public boolean update(Path path, TransitionList transList) throws PrismLangException
	{

		if (valueKnown) {
			return true;
		}

		checkIndx++;
		transLen = path.size();

		if (checkIndx == checkBound) {
			checkIndx = 0;
			updateTmr.start();
			PathFull pf = (PathFull) path;

			Candidate<Pair<State, Integer>> candidate = updateCandidate(pf);

			processTmr.start();
			if (!candidate.equals(lastNonTrivialCnd) && !candidate.isTrivial()) {
				// a new non-trivial candidate is found - increase VN
				// Maxi: visitNo corresponds to k_i in alg 1 in FSMC paper. It can be transformed to a base term and something that is added for every new seen candidate.
				// Maxi: Slightly adapted baseVN for the monitor paper, to account for formula change. Now it is alpha*(i-log epsilon),
				// where alpha = -1/log(1-pmin) and epsilon somehow depends on delta, so I just use delta parameter for it
				// Actually, nothing changed, I think
				cndFound++;
				visitNo = (int) Math.ceil((cndFound - 1.0) * increseVN + baseVN);
				lastNonTrivialCnd = candidate;
				System.out.println("Resets: " + numSamples + " visitNo: " + visitNo + " candSize: " + candidate.states.size()); // + " cand: " + candidate);
			}
			// new exp update: check only up to strength 30. For 30, we check if happy. if happy, terminate exp (hack)
			// if not happy, continue as before, and never reset bold.
			if (cndFound > 0 && candidate.isStrong(Math.min(visitNo,10))) {
				try {
					value = checkProperty(candidate, pf);

					entryToBirth = candidate.getBirthIdx() - candidate.firstentry;
					if(!value){
						//if the value is false, we certainly stop, so we can say valueKnown.
						//however, new exp update: need to check whether bold would really reset now
						if (increseVN>0){//bold continues checking for actual strength
							valueKnown = candidate.isStrong(visitNo);
						}
						else{ //cautious stops
							valueKnown=true;
						}

					}
					else{
						// TODO: If value is true, check whether we certainly are in happy BSCC. Then we can terminate experiment.
						// TODO: But also check this before enourmous bound, or maybe timeout related?

						if(increseVN>0){ //if cautious monitor, this was set to 0, as we don't increase visitnumber. For bold, this has to be greater than 0, so the check makes sense.
							//if bold, set valueKnown to true, as the candidate is strong enough to be certain
							valueKnown = true;
						}
						else{
							//if cautious and the value is true, we can't set valueKnown to true yet, as we aren't certain about the cycle being a BSCC
//							//use the alpha*(i-log epsilon) formula, where alpha = -1/log(1-pmin)
//							int strongness = (int) Math.ceil((-1.0 / ld(simMinProb)) * (cndFound - ld(0.1)));
//							valueKnown = candidate.isStrong(strongness);
							//New exp update: we just use some strength that's good enough to make sure we really are the BSCC
							valueKnown = candidate.isStrong(10);
						}
					}

				} catch (PrismException e) {
					throw new PrismLangException(e.getMessage());
				}

				transLen = ((long) candidate.firstentry);
				candLen = pf.size() - transLen;
			}

			lastIndx = (int) pf.size();
			processTmr.stop();
			updateTmr.stop();
		}

		return valueKnown;
	}

	@Override
	public void updateStats()
	{
		numSamples++;

		canTotal += candLen;
		transTotal += transLen;
		cndFoundTotal += cndFound;
		cndFoundMax = cndFoundMax < cndFound ? cndFound : cndFoundMax;
		entryToBirthTotal += entryToBirth;
		lastRunSteps = candLen+transLen;
		candLen = 0;
		transLen = 0;
		//cndFound = 0;
		entryToBirth = 0;

		if (value) {
			numTrue++;
		}
	}

	@Override
	public Object getCurrentValue()
	{
		return new Boolean(value);
	}

	@Override
	public double getMeanValue()
	{
		return numTrue / (double) numSamples;
	}

	@Override
	public double getVariance()
	{
		// Estimator to the variance (see p.24 of Vincent Nimal's MSc thesis)
		if (numSamples <= 1) {
			return 0.0;
		} else {
			return (numTrue * ((double) numSamples - numTrue) / (numSamples * (numSamples - 1.0)));
		}

		// An alternative, below, would be to use the empirical mean
		// (this is not equivalent (or unbiased) but, asymptotically, is the same)
		//double mean = numTrue / (double) numSamples;
		//return mean * (1.0 - mean);
	}

	@Override
	public double getLikelihoodRatio(double p1, double p0) throws PrismException
	{
		// See Sec 5.3 of Vincent Nimal's MSc thesis for details
		return Math.pow(p1 / p0, numTrue) * Math.pow((1 - p1) / (1 - p0), numSamples - numTrue);
	}

	/**
	 * Check if the candidate satisfies the property.
	 *
	 * @param cnd
	 * @param pf
	 * @return
	 * @throws PrismException
	 */
	private boolean checkProperty(Candidate<Pair<State, Integer>> cnd, PathFull pf) throws PrismException
	{
		constructTmr.start();

		NatBitSet states = NatBitSets.set();
		for (Pair<State, Integer> st : cnd.getStates()) {
			states.set(st.second);
		}

		boolean res = dra.getAcceptance().isBSCCAccepting(states);

		checkTmr.stop();
		return res;
	}

	/**
	 * Set of labels ids for the state
	 *
	 * @param st
	 * @return
	 * @throws PrismLangException
	 */
	private NatBitSet computeLabels(State st) throws PrismLangException
	{
		NatBitSet bs = NatBitSets.set();
		for (int i = 0; i < this.labelExp.length; i++) {
			if (labelExp[i].evaluateBoolean(st)) {
				bs.set(i);
			}
		}

		return bs;
	}

	/**
	 * Returns the candidate for the path.
	 *
	 * @param path
	 * @return
	 * @throws PrismLangException
	 */
	private Candidate<Pair<State, Integer>> updateCandidate(PathFull path) throws PrismLangException
	{
		cndTmr.start();
		Candidate<Pair<State, Integer>> cnd = cndLst.size() > 0 ? cndLst.getLast() : null;

		for (int i = lastIndx + 1; i <= path.size(); i++) {
			State state = path.getState(i);
			NatBitSet labels = computeLabels(state);
			draState = dra.getEdgeDestByLabel(draState, labels);
			Pair<State, Integer> st = new Pair<>(state, draState);

			if (cnd != null && cnd.contains(st)) {
				// state already in the current candidate
				cnd.countState(st);
				cnd.setNonTrivial();
			} else {
				Integer indx = cndMap.get(st);

				if (indx == null) {
					// state seen for the first time
					cnd = new Candidate<>(i, st);
					cndLst.add(cnd);
					cndMap.put(st, cndLst.size() - 1);
				} else {
					// state seen before - create a merge of all candidates from that point
					cnd = new Candidate<>(i, st);

					while (cndLst.size() > indx) {
						Candidate<Pair<State, Integer>> oldCnd = cndLst.removeLast();
						cnd.mergeWith(oldCnd);
					}

					cndLst.addLast(cnd);
					int v = cndLst.size() - 1;

					for (Pair<State, Integer> s : cnd.getStates()) {
						cndMap.put(s, v);
					}

					cnd.setNonTrivial();
				}
			}
		}

		cndTmr.stop();
		return cnd;
	}

	public String toString()
	{
		StringBuilder bldr = new StringBuilder();
		double cndFoundAvg = ((double) cndFoundTotal) / ((double) numSamples);
		bldr.append(this.getClass().getName() + " statistics:\n");
		bldr.append("Visit no. (base, increment):        (" + PrismUtils.formatDouble(5, baseVN) + ", " + PrismUtils.formatDouble(5, increseVN) + ")\n");
		bldr.append("Time for path update:               " + updateTmr + "\n");
		bldr.append("Time for updating candidates:       " + cndTmr + "\n");
		bldr.append("Time for processing candidates:     " + processTmr + "\n");
		bldr.append("Time for DTMC constr.               " + constructTmr + "\n");
		bldr.append("Time for prop. checking.:           " + checkTmr + "\n");
		bldr.append("Number of cnd. (max, avg.):         (" + cndFoundMax + ", " + PrismUtils.formatDouble(4, cndFoundAvg) + ")\n");
		if (countFP) {
			double falseCndPr = ((double) totalFP * 100) / ((double) numSamples);
			bldr.append("No of false positivies:             " + PrismUtils.formatDouble(2, falseCndPr) + "%\n");
		}
		bldr.append("Avg. (transient,BSCC) path len.:    (" + transTotal / numSamples + ", " + canTotal / numSamples + ")\n");
		bldr.append("Avg. entry to birth time:           " + entryToBirthTotal / numSamples + "\n");
		bldr.append("Empirical mean:                     " + getMeanValue() + "\n");

		return bldr.toString();
	}

	public int getNumSamples()
	{
		return numSamples;
	}

	public long getCanTotal()
	{
		return canTotal;
	}

	public long getTransTotal()
	{
		return transTotal;
	}

	public long getLastRunSteps(){
		return lastRunSteps;
	}
}
