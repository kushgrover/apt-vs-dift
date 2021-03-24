/*
 * Copyright (C) 2016  (Salomon Sickert)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package explicit;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.List;

public class ECComputerFast
{
	// TODO Merge / replace ECComputer

	public static class DecompositionResult {
		public final List<NatBitSet> sccs;
		public final List<NatBitSet> bsccs;

		public DecompositionResult(List<NatBitSet> sccs, List<NatBitSet> bsccs)
		{
			this.sccs = sccs;
			this.bsccs = bsccs;
		}
	}

	public static DecompositionResult computeSCCs(Model model)
	{
		NatBitSet states = NatBitSets.boundedFilledSet(model.getNumStates());
		return computeSCCs(model, states);
	}

	public static DecompositionResult computeSCCs(Model model, NatBitSet states)
	{
		return Tarjan.getSCCs(new Tarjan.DefaultTarjanModel(model), states);
	}

	public static List<MEC> computeMECs(NondetModel model)
	{
		return computeMECs(model, null);
	}

	/**
	 * Computes maximal end components of a nondeterministic model such as an MDP.
	 * The implementation uses a variation of the algorithm from p.48 of:
	 * Luca de Alfaro. Formal Verification of Probabilistic Systems. Ph.D. thesis, Stanford University (1997)
	 *
	 * @param model   The model
	 * @param preMECs if null all reachable states are considered.
	 * @return A list of all MECs obtained from refining the preMECs.
	 */
	public static List<MEC> computeMECs(NondetModel model, Collection<MEC> preMECs)
	{
		List<MEC> mecs = new ArrayList<>();
		Deque<MEC> workList = new ArrayDeque<>();

		if (preMECs == null) {
			NatBitSet bs = NatBitSets.boundedFilledSet(model.getNumStates());
			workList.add(MEC.createMEC(model, bs));
		} else {
			workList.addAll(preMECs);
		}

		List<MEC> refinedMECs = new ArrayList<>();
		Tarjan.SCCConsumer preMecConsumer = preMecStates -> {
			if (preMecStates.isEmpty()) {
				return true;
			}
			MEC preMec = MEC.createMEC(model, preMecStates);
			if (preMec.states.isEmpty()) {
				return true;
			}
			assert !refinedMECs.contains(preMec);
			refinedMECs.add(preMec);
			return true;
		};

		while (!workList.isEmpty()) {
			MEC mec = workList.remove();

			refinedMECs.clear();
			Tarjan.runTarjan(new Tarjan.NondetMecTarjanModel(model, mec), preMecConsumer, mec.states);
			if (!refinedMECs.isEmpty()) {
				if (refinedMECs.size() == 1) {
					MEC refinedMEC = refinedMECs.get(0);
					if (mec.equals(refinedMEC)) {
						assert !mecs.contains(refinedMEC);
						mecs.add(refinedMEC);
						continue;
					}
				}
				workList.addAll(refinedMECs);
			}
		}

		return mecs;
	}

	public static List<NatBitSet> computeBSCCs(Model model, NatBitSet states)
	{
		if (states == null) {
			int numStates = model.getNumStates();
			states = NatBitSets.boundedFilledSet(numStates, numStates);
		}
		return Tarjan.getBSCCs(new Tarjan.DefaultTarjanModel(model), states);
	}

	public static List<NatBitSet> computeBSCCs(Model model)
	{
		return computeBSCCs(model, NatBitSets.boundedFilledSet(model.getNumStates()));
	}
}
