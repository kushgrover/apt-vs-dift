package simulator.sampler;

import cern.colt.Timer;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import parser.State;
import parser.Values;
import parser.VarList;
import parser.ast.Expression;
import parser.ast.ExpressionTemporal;
import parser.ast.ModulesFile;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismUtils;
import simulator.Choice;
import simulator.Path;
import simulator.PathFull;
import simulator.TransitionList;
import simulator.UnboundedSimulationParameters;
import simulator.Updater;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Deque;
import java.util.HashSet;
import java.util.Set;

/**
 * Sampler checks if the trace satisfies a simple until property by identifying a BSCC candidate in the trace.
 *
 * @author Przemysław Daca
 */
public class SamplerCandidateUntil extends SamplerBoolean
{
	// the initial visit number
	private final double baseVN;
	// whether we can looking into the model or only sample
	private final boolean whitebox;
	private Updater updater;
	// the max. probability the sampler identifies a false candidate
	private double simFalseCnd;
	// the minimal visit number (VN) of status in the candidate
	private int visitNo;
	// how many times we increase VN after a new non-trivial candidate is found
	private double increseVN;
	// how often to check for a VN-candidate
	private int checkBound;
	private int checkIndx;

	// CHECKING CANDIDATES
	// last inspected index	(-1 if none)
	private int lastIndx;
	//private Candidate candidate;
	private Candidate<State> lastNonTrivialCnd;
	// sequence of non-subsummed candidates  on the path
	private Deque<Candidate<State>> cndLst;
	// state -> index in cndLst
	private Object2IntMap<State> cndMap;

	// STATISTICS
	// count the number of candidates that were false positives
	private boolean countFP = true;

	private Timer updateTmr = new Timer();
	private Timer cndTmr = new Timer();
	private Timer processTmr = new Timer();

	private long candLen;
	private long transLen;
	private long canTotal;
	private long transTotal;
	private int totalFP;
	private int entryToBirth;
	private int cndFound;
	private int cndFoundMax;
	private int cndFoundTotal;
	private Expression left;
	private Expression right;
	private int entryToBirthTotal;
	private int candidateDiscovered;
	private int largestCnd;

	public SamplerCandidateUntil(ModulesFile mf, ExpressionTemporal expr, Values constantValues, UnboundedSimulationParameters usp) throws PrismException
	{
		this.simFalseCnd = usp.getSimFalseCnd();
		double ratio = usp.getSimRatio();
		this.checkBound = usp.getSimCheckBound();
		this.whitebox = usp.isSimWhitebox();

		this.cndMap = new Object2IntOpenHashMap<>();
		this.cndLst = new ArrayDeque<>();
		this.baseVN = Math.log((1.0 - ratio) * this.simFalseCnd) / Math.log(1 - usp.getSimMinProb()) + 1;
		this.increseVN = Math.log(ratio) / Math.log(1 - usp.getSimMinProb());
		this.visitNo = (int) Math.ceil(baseVN);

		if (expr.getOperator() != ExpressionTemporal.P_U) {
			throw new PrismException("Error creating Sampler");
		}
		left = expr.getOperand1();
		right = expr.getOperand2();

		mf = (ModulesFile) mf.deepCopy().replaceConstants(constantValues).simplify();
		PrismComponent pc = new PrismComponent();
		VarList varList = mf.createVarList();
		this.updater = new Updater(mf, varList, pc);
	}

	@Override
	public void reset()
	{
		super.reset();

		visitNo = (int) Math.ceil(baseVN);

		checkIndx = 0;
		lastIndx = -1;
		lastNonTrivialCnd = null;
		cndMap.clear();
		cndLst.clear();
	}

	@Override
	public void resetStats()
	{
		super.resetStats();
		candLen = 0;
		transLen = 0;
		cndFound = 0;
	}

	@Override
	public boolean update(Path path, TransitionList transList) throws PrismLangException
	{

		if (valueKnown) {
			return true;
		}

		updateTmr.start();

		checkIndx++;
		transLen = path.size();

		if (checkIndx == checkBound) {
			checkIndx = 0;
			PathFull pf = (PathFull) path;

			Candidate<State> candidate = updateCandidate(pf);

			largestCnd = Math.max(largestCnd, candidate.getStates().size());

			processTmr.start();
			if (!candidate.equals(lastNonTrivialCnd) && !candidate.isTrivial()) {
				// a new non-trivial candidate is found - increase VN
				cndFound++;
				visitNo = (int) Math.ceil((cndFound - 1.0) * increseVN + baseVN);
				lastNonTrivialCnd = candidate;
			}

			if (whitebox && !candidate.isTrivial() && !candidate.isIncomplete()) {
				try {
					if (verifyCandidate(candidate)) {
						// we hit a candidate, without reaching the target state
						valueKnown = true;
						value = false;
					} else {
						totalFP++;
						candidate.setIncomplete();
					}
				} catch (PrismException e) {
					throw new PrismLangException(e.getMessage());
				}

				transLen = ((long) candidate.firstentry);
				candLen = pf.size() - transLen;
			}

			if (!whitebox && cndFound > 0 && candidate.isStrong(visitNo)) {
				try {
					if (countFP && !verifyCandidate(candidate)) {
						totalFP++;
					}

					// we hit a candidate, without reaching the target state
					valueKnown = true;
					value = false;
					entryToBirth = candidate.getBirthIdx() - candidate.firstentry;
					candidateDiscovered++;
				} catch (PrismException e) {
					throw new PrismLangException(e.getMessage());
				}

				transLen = ((long) candidate.getBirthIdx());
				candLen = pf.size() - transLen;
			}

			lastIndx = (int) pf.size();
			processTmr.stop();
		}

		//  for simple until properties, we can simply check if the LHS is satisfied or RHS is violated
		if (!valueKnown) {
			State state = path.getCurrentState();
			if (right.evaluateBoolean(state)) {
				valueKnown = true;
				value = true;
			} else if (!left.evaluateBoolean(state)) {
				valueKnown = true;
				value = false;
			}
		}
		updateTmr.stop();

		return valueKnown;
	}

	/**
	 * Check if the candidate is really a BSCC.
	 *
	 * @param cand
	 * @throws PrismException
	 */
	private boolean verifyCandidate(Candidate<State> cand) throws PrismException
	{
		// check if the candidate is closed under transition relation
		Set<State> states = cand.getStates();
		for (State st : states) {
			Collection<State> succs = getSucessors(st);

			if (!states.containsAll(succs)) {
				return false;
			}
		}
		return true;
	}

	private Collection<State> getSucessors(State st) throws PrismException
	{
		TransitionList tlist = new TransitionList();
		updater.calculateTransitions(st, tlist);

		Set<State> succs = new HashSet<>();

		for (int i = 0; i < tlist.getNumChoices(); i++) {
			Choice choice = tlist.getChoice(i);

			for (int j = 0; j < choice.size(); j++) {
				State succ = choice.computeTarget(j, st);
				succs.add(succ);
			}
		}

		return succs;
	}

	@Override
	public void updateStats()
	{
		super.updateStats();

		canTotal += candLen;
		transTotal += transLen;
		cndFoundTotal += cndFound;
		cndFoundMax = cndFoundMax < cndFound ? cndFound : cndFoundMax;
		entryToBirthTotal += entryToBirth;
		candLen = 0;
		transLen = 0;
		cndFound = 0;
	}

	/**
	 * Returns the candidate for the path.
	 *
	 * @param path
	 * @return
	 */
	private Candidate<State> updateCandidate(PathFull path)
	{
		cndTmr.start();
		Candidate<State> cnd = cndLst.isEmpty() ? null : cndLst.getLast();

		for (int i = lastIndx + 1; i <= path.size(); i++) {
			State st = path.getState(i);
			if (cnd != null && cnd.contains(st)) {
				// state already in the current candidate
				if (!whitebox) {
					cnd.countState(st);
				}
				cnd.setNonTrivial();
			} else {
				int indx = cndMap.getOrDefault(st, -1);
				if (indx == -1) {
					// state seen for the first time
					cnd = new Candidate<>(i, st);
					cndLst.add(cnd);
					cndMap.put(st, cndLst.size() - 1);
				} else {
					// state seen before - create a merge of all candidates from that point
					cnd = new Candidate<>(i, st);

					while (cndLst.size() > indx) {
						Candidate<State> oldCnd = cndLst.removeLast();
						cnd.mergeWith(oldCnd);
					}

					cndLst.addLast(cnd);
					int v = cndLst.size() - 1;

					for (State state : cnd.getStates()) {
						cndMap.put(state, v);
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
		bldr.append(this.getClass().getName()).append(" statistics:\n");
		bldr.append("Visit no. (base, increment):        (").append(PrismUtils.formatDouble(5, baseVN)).append(", ")
				.append(PrismUtils.formatDouble(5, increseVN)).append(")\n");
		bldr.append("Time for path update:               ").append(updateTmr).append("\n");
		bldr.append("Time for updating candidates:       ").append(cndTmr).append("\n");
		bldr.append("Time for processing candidates:     ").append(processTmr).append("\n");
		bldr.append("Number of cnd. (max, avg.):         (").append(cndFoundMax).append(", ").append(PrismUtils.formatDouble(4, cndFoundAvg)).append(")\n");
		if (!whitebox && countFP) {
			double falseCndPr = ((double) totalFP * 100.0) / ((double) numSamples);
			bldr.append("No of false positivies:             ").append(PrismUtils.formatDouble(2, falseCndPr)).append("%\n");
		}
		if (whitebox) {
			bldr.append("No of false candidates:             ").append(totalFP).append("\n");
		}
		bldr.append("Size of largest candidate:          ").append(largestCnd).append("\n");
		bldr.append("Avg. (transient,BSCC) path len.:    (").append(((double) transTotal) / ((double) numSamples)).append(", ")
				.append(((double) canTotal) / ((double) numSamples)).append(")\n");
		if (candidateDiscovered > 0) {
			bldr.append("Avg. entry to birth time:           ").append(entryToBirthTotal / candidateDiscovered).append("\n");
		}
		bldr.append("Empirical mean:                     ").append(getMeanValue()).append("\n");


		return bldr.toString();
	}
}
