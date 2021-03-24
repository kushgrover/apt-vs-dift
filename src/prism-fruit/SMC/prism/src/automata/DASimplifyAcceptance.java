package automata;

import de.tum.in.naturals.set.NatBitSet;

import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.IntIterator;
import prism.PrismComponent;
import prism.PrismException;
import explicit.SCCComputer;
import acceptance.AcceptanceOmega;
import acceptance.AcceptanceRabin;
import acceptance.AcceptanceReach;
import acceptance.AcceptanceType;
import acceptance.AcceptanceRabin.RabinPair;

public class DASimplifyAcceptance
{

	/**
	 * Tries to simplify the acceptance condition of the deterministic automaton.
	 * Note that the passed parameter {@code da} may be destroyed by this function.
	 * @param parent the calling PrismComponent (for SCC computer)
	 * @param da the DA to be simplified (may be destroyed)
	 * @param allowedAcceptance the allowed acceptance types
	 */
	@SuppressWarnings("unchecked")
	public static DA<NatBitSet, ? extends AcceptanceOmega> simplifyAcceptance(PrismComponent parent, DA<NatBitSet, ? extends AcceptanceOmega> da, AcceptanceType... allowedAcceptance)
			throws PrismException
	{
		// Simplifications for DRAs
		if (da.getAcceptance() instanceof AcceptanceRabin) {
			DA<NatBitSet, AcceptanceRabin> dra = (DA<NatBitSet, AcceptanceRabin>) da;
			// K_i states that do not occur in a (non-trivial) SCC of the DRA may as well be removed
			SCCComputer sccComp = explicit.SCCComputer.createSCCComputer(parent, new LTSFromDA(da));
			sccComp.computeBSCCs();
			NatBitSet trivial = sccComp.getNotInSCCs();
			for (RabinPair pair : dra.getAcceptance()) {
				if (pair.getK().intersects(trivial)) {
					pair.getK().andNot(trivial);
				}
			}
			// See if the DRA is actually a DFA
			if (AcceptanceType.contains(allowedAcceptance, AcceptanceType.REACH) && isDFA(dra)) {
				// we can switch to AcceptanceReach
				AcceptanceReach reachAcceptance = new AcceptanceReach(getDFAGoalStatesForRabin(dra.getAcceptance()));
				DA.switchAcceptance(dra, reachAcceptance);
				da = dra;
			}
		}
		return da;
	}

	/**
	 * Is this Rabin automaton actually a finite automaton? This check is done syntactically:
	 * it returns true if every transition from a K_i state goes to another K_i state.
	 * We also require that there are no L_i states overlapping with any K_j states.
	 */
	public static boolean isDFA(DA<NatBitSet, AcceptanceRabin> dra)
	{
		AcceptanceRabin acceptance = dra.getAcceptance();
		// Compute potential set of goal states as the union of all K_i sets
		NatBitSet goalStates = getDFAGoalStatesForRabin(acceptance);

		// Make sure there are no L_i states in the goal states for any i
		for (RabinPair anAcceptance : acceptance) {
			if (goalStates.intersects(anAcceptance.getL()))
				return false;
		}
		// Check if every transition from a goal state goes to another goal state
		IntIterator iterator = goalStates.iterator();
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			int m = dra.getNumEdges(i);
			for (int j = 0; j < m; j++) {
				if (!goalStates.contains(dra.getEdgeDest(i, j)))
					return false;
			}
		}
		return true;
	}

	/**
	 * Get the union of the K_i states of a Rabin acceptance condition.
	 */
	public static NatBitSet getDFAGoalStatesForRabin(AcceptanceRabin acceptance)
	{
		// Compute set of goal states as the union of all K_i sets
		NatBitSet goalStates = NatBitSets.set();
		int n = acceptance.size();
		for (RabinPair anAcceptance : acceptance) {
			goalStates.or(anAcceptance.getK());
		}
		return goalStates;
	}
}
