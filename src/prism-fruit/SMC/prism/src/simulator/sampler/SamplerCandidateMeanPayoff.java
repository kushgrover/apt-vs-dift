package simulator.sampler;

import cern.colt.Timer;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import explicit.DTMCModelChecker;
import explicit.DTMCSimple;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import parser.State;
import parser.Values;
import parser.VarList;
import parser.ast.ExpressionReward;
import parser.ast.ModulesFile;
import parser.ast.RewardStruct;
import prism.Pair;
import prism.PrismComponent;
import prism.PrismDevNullLog;
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
 * Sampler computes a mean-payoff value by identifying a BSCC candidate in the trace.
 *
 * @author Przemysław Daca
 */
public class SamplerCandidateMeanPayoff extends SamplerDouble
{
	private final Values constantValues;
	private final Updater updater;
	// the initial visit number
	private final double baseVN;
	// rewards
	private final int rewardStructIndex;
	private final int numRS;
	private final Timer updateTmr = new Timer();
	private final Timer cndTmr = new Timer();
	private final Timer processTmr = new Timer();
	private ModulesFile mf;
	// the minimal visit number (VN) of status in the candidate
	private int visitNo;
	private double visitNoEst;
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
	private Integer cndEdges;
	// STATISTICS
	// count the number of candidates that were false positives
	private boolean countFP = true;
	private long candLen;
	private long transLen;
	private int totalFP;
	private int cndFound;
	private int cndFoundTotal;
	private int transTotal;
	private int cndFoundMax;
	private int canTotal;

	// TODO do sth with the magic values
	private double width_div = 0.08;
	private double div_probability = 0.01;
	private double minprob;
	private DTMCModelChecker mcDTMC;

	public SamplerCandidateMeanPayoff(ModulesFile mf, Values constantValues, int rewardStructIndex, ExpressionReward expr, UnboundedSimulationParameters usp)
			throws PrismException
	{
		this.mf = mf;
		this.constantValues = constantValues;
		double simFalseCnd = usp.getSimFalseCnd();
		double ratio = usp.getSimRatio();
		this.checkBound = usp.getSimCheckBound();
		this.cndMap = new Object2IntOpenHashMap<>();
		this.cndLst = new ArrayDeque<>();
		this.minprob = usp.getSimMinProb();
		this.baseVN = Math.log((1.0 - ratio) * simFalseCnd) / Math.log(1.0 - this.minprob) + 1;
		this.increseVN = Math.log(ratio) / Math.log(1.0 - this.minprob);
		visitNo = (int) Math.ceil(baseVN);

		this.mcDTMC = new DTMCModelChecker(null);
		this.mcDTMC.setLog(new PrismDevNullLog());
		this.rewardStructIndex = rewardStructIndex;
		this.numRS = mf.getNumRewardStructs();

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
		visitNoEst = 0;

		checkIndx = 0;
		lastIndx = -1;
		lastNonTrivialCnd = null;
		cndEdges = null;
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
	public void updateStats()
	{
		super.updateStats();

		canTotal += candLen;
		transTotal += transLen;
		cndFoundTotal += cndFound;
		cndFoundMax = cndFoundMax < cndFound ? cndFound : cndFoundMax;
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

			processTmr.start();
			if (!candidate.equals(lastNonTrivialCnd) && !candidate.isTrivial()) {
				// a new non-trivial candidate is found - increase VN
				cndFound++;
				visitNo = (int) Math.ceil((cndFound - 1.0) * increseVN + baseVN);
				lastNonTrivialCnd = candidate;
				cndEdges = null;
			}

			if (cndFound > 0 && candidate.isStrong(visitNo)) {

				if (cndEdges == null) {
					cndEdges = countEdges(candidate, pf);

					int size = candidate.getStates().size();

					if (size == 1) {
						visitNoEst = 0;
					} else {
						double delta_pr = Math.pow(width_div + 1.0, 1.0 / (2.0 * size)) - 1.0;
						delta_pr = delta_pr * this.minprob;
						visitNoEst = Math.log((2.0 * (cndEdges - size)) / div_probability) / Math.log(Math.E);
						visitNoEst = visitNoEst / (2 * delta_pr * delta_pr);
					}

					if (visitNoEst > Integer.MAX_VALUE) {
						throw new PrismLangException("Sample too long");
					}

				}

				if (candidate.isStrong(visitNo + ((int) visitNoEst))) {
					try {
						if (countFP && !verifyCandidate(candidate)) {
							totalFP++;
						}

						valueKnown = true;
						this.value = computeMP(candidate, pf);
					} catch (PrismException e) {
						throw new PrismLangException(e.getMessage());
					}
				}

				transLen = candidate.getBirthIdx();
				candLen = pf.size() - transLen;
			}

			lastIndx = (int) pf.size();
			processTmr.stop();
		}

		//  for simple until properties, we can simply check if the LHS is satisfied or RHS is violated
		updateTmr.stop();

		return valueKnown;
	}

	private Integer countEdges(Candidate<State> candidate, PathFull pf)
	{
		Set<Pair<State, State>> edgesSet = new HashSet<>();
		int i = candidate.getBirthIdx();
		State prev = pf.getState(i);
		i++;

		while (((long) i) <= pf.size()) {
			State st = pf.getState(i);
			Pair<State, State> pair = new Pair<>(prev, st);
			edgesSet.add(pair);
			i++;
			prev = st;
		}

		return edgesSet.size();
	}

	private double computeMP(Candidate<State> candidate, PathFull pf) throws PrismException
	{
		int size = candidate.getStates().size();
		Object2IntMap<State> stMap = new Object2IntOpenHashMap<>();
		State[] stArr = new State[size];
		double[] rewards = new double[size];

		int i = 0;
		for (State st : candidate.getStates()) {
			stMap.put(st, i);
			stArr[i] = st;
			double[] store = new double[numRS];
			calculateStateRewards(st, store);
			rewards[i] = store[this.rewardStructIndex];
			i++;
		}

		// estimate probabilities
		long[][] probs = new long[size][size];

		int prevIdx, idx;
		i = candidate.getBirthIdx();
		State st = pf.getState(i);
		idx = stMap.getInt(st);

		while (i < pf.size()) {
			prevIdx = idx;
			i++;
			st = pf.getState(i);
			idx = stMap.getInt(st);
			probs[prevIdx][idx]++;
		}

		// build model
		NatBitSet bs = NatBitSets.set();
		DTMCSimple dtmc = new DTMCSimple();
		for (i = 0; i < size; i++) {
			dtmc.addState();
			bs.set(i);
		}

		for (i = 0; i < size; i++) {
			st = stArr[i];
			int cnt = candidate.getCount(st);

			for (int j = 0; j < size; j++) {
				double pr = ((double) probs[i][j]) / cnt;
				if (pr > 0) {
					dtmc.addToProbability(i, j, pr);
				}
			}
		}
		dtmc.addInitialState(0);

		double ss[] = new double[size];
		mcDTMC.computeSteadyStateProbsForBSCC(dtmc, bs, ss);

		// finally, compute the mean-pay off
		double mp = 0;
		for (i = 0; i < size; i++) {
			mp += ss[i] * rewards[i];
		}

		return mp;
	}

	private void calculateStateRewards(State state, double[] store) throws PrismLangException
	{
		int i, j, n;
		double d;
		RewardStruct rw;
		for (i = 0; i < this.numRS; i++) {
			rw = this.mf.getRewardStruct(i);
			n = rw.getNumItems();
			d = 0.0;
			for (j = 0; j < n; j++) {
				if (!rw.getRewardStructItem(j).isTransitionReward()) {
					if (rw.getStates(j).evaluateBoolean(state)) {
						d += rw.getReward(j).evaluateDouble(this.constantValues, state);
					}
				}
			}
			store[i] = d;
		}
	}

	/**
	 * Check if the candidate is really a BSCC.
	 *
	 * @param cand
	 * @return
	 * @throws PrismException
	 */
	private boolean verifyCandidate(Candidate<State> cand) throws PrismException
	{
		// check if the candidate is closed under transition relatio		bldr.append("Time for prop. checking.:           "+checkTmr+"\n");n
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

	/**
	 * Returns the candidate for the path.
	 *
	 * @param path
	 * @return
	 */
	private Candidate<State> updateCandidate(PathFull path)
	{
		cndTmr.start();
		Candidate<State> cnd = !cndLst.isEmpty() ? cndLst.getLast() : null;

		for (int i = lastIndx + 1; i <= path.size(); i++) {
			State st = path.getState(i);

			if (cnd != null && cnd.contains(st)) {
				// state already in the current candidate
				cnd.countState(st);
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
		if (countFP) {
			double falseCndPr = ((double) totalFP * 100) / ((double) numSamples);
			bldr.append("No of false positivies:             ").append(PrismUtils.formatDouble(2, falseCndPr)).append("%\n");
		}
		bldr.append("Avg. (transient,BSCC) path len.:    (").append(transTotal / numSamples).append(", ").append(canTotal / numSamples).append(")\n");
		bldr.append("Empirical mean:                     ").append(getMeanValue()).append("\n");

		return bldr.toString();
	}

}
