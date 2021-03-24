package explicit;

import de.tum.in.naturals.set.BoundedNatBitSet;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntIterator;

import java.util.Iterator;
import java.util.List;
import java.util.function.IntConsumer;

public class DTMCCollapsed extends DTMCSimple
{
	final BoundedNatBitSet bsccStates;
	final Int2IntMap originalToCollapsedStateMap;

	/**
	 * Collapses BSCCs into single states.
	 *
	 * @param dtmc  DTMC to be collapsed
	 * @param bsccs Precomputed BSCCs
	 */
	public DTMCCollapsed(DTMC dtmc, NatBitSet subset, List<NatBitSet> bsccs)
	{
		final int numStates = subset == null ? dtmc.getNumStates() : subset.size();

		// Map from every state to its representative in the new DTMC
		// All states belonging to one BSCC will be mapped to a single representative
		// All transient states will have their own unique representative
		originalToCollapsedStateMap = new Int2IntArrayMap(numStates);

		// Compute bscc states
		bsccStates = NatBitSets.boundedSet(dtmc.getNumStates());
		bsccs.forEach(bsccStates::or);

		// Compute transient states
		NatBitSet transientStates;
		if (subset == null) {
			transientStates = bsccStates.complement();
		} else {
			transientStates = subset.clone();
			transientStates.andNot(bsccStates);
		}

		int i = 0;
		IntIterator transientStateIterator = transientStates.iterator();
		while (transientStateIterator.hasNext()) {
			int transientState = transientStateIterator.nextInt();
			originalToCollapsedStateMap.put(transientState, i);
			i++;
		}

		final int numTransientStates = transientStates.size();
		final int numBsccs = bsccs.size();

		// Add transient states, BSCC representatives and the special + and - states
		this.numStates = numTransientStates + numBsccs;

		super.initialise(this.numStates);
		numTransitions = 0;

		// Create new states B_x for each BSCC B
		for (int bsccIndex = 0; bsccIndex < numBsccs; bsccIndex++) {
			final NatBitSet bscc = bsccs.get(bsccIndex);
			final int currentBsccRepresentativeIndex = numTransientStates + bsccIndex;

			// For every state in the BSCC, map it to its representative
			bscc.forEach((IntConsumer) s -> originalToCollapsedStateMap.put(s, currentBsccRepresentativeIndex));
		}

		// For each transient state, if it has a transition going into some state of B, then change its destination to B_x
		int collapsedState = 0;
		transientStateIterator = transientStates.iterator();
		while (transientStateIterator.hasNext()) {
			int transientState = transientStateIterator.nextInt();

			Distribution collapsedDistribution = new Distribution();
			Iterator<Int2DoubleMap.Entry> transitionIterator = dtmc.getTransitionsIterator(transientState);
			// For every transition from that state (target, probability),
			// If deadlocked state
			if (transitionIterator == null) {
				addDeadlockState(transientState);
				continue;
			}

			// otherwise add (collapsed representative of target, probability) to distribution
			transitionIterator.forEachRemaining(transition ->
					collapsedDistribution.add(originalToCollapsedStateMap.get(transition.getIntKey()), transition.getDoubleValue()));
			trans.set(collapsedState, collapsedDistribution);
			collapsedState++;
			numTransitions++;
		}
	}

	public int getCollapsedStateIndex(int originalState)
	{
		return originalToCollapsedStateMap.get(originalState);
	}
}
