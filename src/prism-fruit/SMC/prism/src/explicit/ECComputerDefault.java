//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	* Mateusz Ujma <mateusz.ujma@cs.ox.ac.uk> (University of Oxford)
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

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import prism.PrismComponent;
import prism.PrismException;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.IntConsumer;

/**
 * Explicit maximal end component computer for a nondeterministic model such as an MDP.
 * Implements the algorithm from p.48 of:
 * Luca de Alfaro. Formal Verification of Probabilistic Systems. Ph.D. thesis, Stanford University (1997)
 */
public class ECComputerDefault extends ECComputer
{
	/** The model to compute (M)ECs for **/
	private final NondetModel model;

	/** Computed list of MECs **/
	private List<NatBitSet> mecs = new ArrayList<>();

	/**
	 * Build (M)EC computer for a given model.
	 */
	public ECComputerDefault(PrismComponent parent, NondetModel model) throws PrismException
	{
		super(parent);
		this.model = model;
	}

	/**
	 * Build (M)EC computer for a given model.
	 */
	public ECComputerDefault(NondetModel model) throws PrismException
	{
		super();
		this.model = model;
	}

	// Methods for ECComputer interface

	@Override
	public void computeMECStates() throws PrismException
	{
		mecs = findEndComponents(null, null);
	}

	@Override
	public void computeMECStates(NatBitSet restrict) throws PrismException
	{
		mecs = findEndComponents(restrict, null);
	}

	@Override
	public void computeMECStates(NatBitSet restrict, NatBitSet accept) throws PrismException
	{
		mecs = findEndComponents(restrict, accept);
	}

	@Override
	public List<NatBitSet> getMECStates()
	{
		return mecs;
	}

	// Computation

	/**
	 * Find all accepting maximal end components (MECs) in the submodel obtained
	 * by restricting this one to the set of states {@code restrict},
	 * where acceptance is defined as those which intersect with {@code accept}.
	 * If {@code restrict} is null, we look at the whole model, not a submodel.
	 * If {@code accept} is null, the acceptance condition is trivially satisfied.
	 * @param restrict NatBitSet for the set of states to restrict to
	 * @param accept NatBitSet for the set of accepting states
	 * @return a list of NatBitSets representing the MECs
	 */
	private List<NatBitSet> findEndComponents(NatBitSet restrict, NatBitSet accept) throws PrismException
	{
		// If restrict is null, look within set of all reachable states
		if (restrict == null) {
			restrict = NatBitSets.boundedFilledSet(model.getNumStates());
		}
		// Initialise L with set of all states to look in (if non-empty)
		List<NatBitSet> L = new ArrayList<>();
		if (restrict.isEmpty())
			return L;
		L.add(restrict);
		// Find MECs
		boolean changed = true;
		while (changed) {
			NatBitSet E = L.remove(0);
			SubNondetModel submodel = restrict(model, E);
			List<NatBitSet> sccs = translateStates(submodel, computeSCCs(submodel));
			L = replaceEWithSCCs(L, E, sccs);
			changed = canLBeChanged(L, E);
		}
		// Filter and return those that contain a state in accept
		if (accept != null) {
			int i = 0;
			while (i < L.size()) {
				if (!L.get(i).intersects(accept)) {
					L.remove(i);
				} else {
					i++;
				}
			}
		}
		return L;
	}

	private final Set<NatBitSet> processedSCCs = new HashSet<>();

	private boolean canLBeChanged(List<NatBitSet> L, NatBitSet E)
	{
		processedSCCs.add(E);
		for (NatBitSet aL : L) {
			if (!processedSCCs.contains(aL)) {
				return true;
			}
		}
		return false;
	}

	private List<NatBitSet> replaceEWithSCCs(List<NatBitSet> L, NatBitSet E, List<NatBitSet> sccs)
	{
		if (!sccs.isEmpty()) {
			List<NatBitSet> toAdd = new ArrayList<>();
			for (NatBitSet scc : sccs) {
				if (!L.contains(scc)) {
					toAdd.add(scc);
				}
			}
			if (!toAdd.isEmpty()) {
				L.addAll(toAdd);
			}
		}
		return L;
	}

	private SubNondetModel restrict(NondetModel model, NatBitSet states)
	{
		Int2ObjectMap<NatBitSet> actions = new Int2ObjectOpenHashMap<>();
		NatBitSet initialStates = NatBitSets.boundedSet(model.getNumStates());
		initialStates.set(states.firstInt());

		boolean changed = true;
		while (changed) {
			changed = false;
			actions.clear();
			for (int i = 0; i < model.getNumStates(); i++) {
				NatBitSet act = NatBitSets.set();
				if (states.contains(i)) {
					for (int j = 0; j < model.getNumChoices(i); j++) {
						if (model.allSuccessorsInSet(i, j, states)) {
							act.set(j);
						}
					}
					if (act.isEmpty()) {
						states.clear(i);
						changed = true;
					}
					actions.put(i, act);
				}
			}
		}

		return new SubNondetModel(model, states, actions, initialStates);
	}

	private List<NatBitSet> computeSCCs(NondetModel model) throws PrismException
	{
		SCCComputer sccc = SCCComputer.createSCCComputer(this, model);
		sccc.computeSCCs();
		return sccc.getSCCs();
	}

	private List<NatBitSet> translateStates(SubNondetModel model, List<NatBitSet> sccs)
	{
		List<NatBitSet> r = new ArrayList<>();
		for (NatBitSet set : sccs) {
			NatBitSet set2 = NatBitSets.set();
			r.add(set2);
			set.forEach((IntConsumer) j -> set2.set(model.translateState(j)));
		}
		return r;
	}

	private boolean isMEC(NatBitSet b)
	{
		if (b.isEmpty()) {
			return false;
		}

		IntIterator iterator = b.iterator();
		while (iterator.hasNext()) {
			int state = iterator.nextInt();
			boolean atLeastOneAction = false;
			for (int i = 0; i < model.getNumChoices(state); i++) {
				if (model.allSuccessorsInSet(state, i, b)) {
					atLeastOneAction = true;
				}
			}
			if (!atLeastOneAction) {
				return false;
			}
		}

		return true;
	}
}
