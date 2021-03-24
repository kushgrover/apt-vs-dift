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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import prism.PrismComponent;
import prism.PrismException;

/**
 * Explicit maximal end component computer for a nondeterministic model such as an MDP.
 * Implements the algorithm from p.48 of:
 * Luca de Alfaro. Formal Verification of Probabilistic Systems. Ph.D. thesis, Stanford University (1997)
 */
public class ECComputerDefault extends ECComputer
{
	/** The model to compute (M)ECs for **/
	private NondetModel model;

	/** Computed list of MECs **/
	private List<BitSet> mecs = new ArrayList<BitSet>();

	/**
	 * Build (M)EC computer for a given model.
	 */
	public ECComputerDefault(PrismComponent parent, NondetModel model) throws PrismException
	{
		super(parent);
		this.model = model;
	}

	// Methods for ECComputer interface

	@Override
	public void computeMECStates() throws PrismException
	{
		mecs = findEndComponents(null, null);
	}

	@Override
	public void computeMECStates(BitSet restrict) throws PrismException
	{
		mecs = findEndComponents(restrict, null);
	}

	@Override
	public void computeMECStates(BitSet restrict, BitSet accept) throws PrismException
	{
		mecs = findEndComponents(restrict, accept);
	}

	@Override
	public List<BitSet> getMECStates()
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
	 * @param restrict BitSet for the set of states to restrict to
	 * @param accept BitSet for the set of accepting states
	 * @return a list of BitSets representing the MECs
	 */
	private List<BitSet> findEndComponents(BitSet restrict, BitSet accept) throws PrismException
	{
		// If restrict is null, look within set of all reachable states
		if (restrict == null) {
			restrict = new BitSet();
			restrict.set(0, model.getNumStates());
		}
		// Initialise L with set of all states to look in (if non-empty)
		List<BitSet> L = new ArrayList<BitSet>();
		if (restrict.isEmpty())
			return L;
		L.add(restrict);
		// Find MECs
		boolean changed = true;
		while (changed) {
			changed = false;
			BitSet E = L.remove(0);
			SubNondetModel submodel = restrict(model, E);
			List<BitSet> sccs = translateStates(submodel, computeSCCs(submodel));
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

	private Set<BitSet> processedSCCs = new HashSet<BitSet>();

	private boolean canLBeChanged(List<BitSet> L, BitSet E)
	{
		processedSCCs.add(E);
		for (int i = 0; i < L.size(); i++) {
			if (!processedSCCs.contains(L.get(i))) {
				return true;
			}
		}
		return false;
	}

	private List<BitSet> replaceEWithSCCs(List<BitSet> L, BitSet E, List<BitSet> sccs)
	{
		if (sccs.size() > 0) {
			List<BitSet> toAdd = new ArrayList<BitSet>();
			for (int i = 0; i < sccs.size(); i++) {
				if (!L.contains(sccs.get(i))) {
					toAdd.add(sccs.get(i));
				}
			}
			if (toAdd.size() > 0) {
				L.addAll(toAdd);
			}
		}
		return L;
	}

	private SubNondetModel restrict(NondetModel model, BitSet states)
	{
		Map<Integer, BitSet> actions = new HashMap<Integer, BitSet>();
		BitSet initialStates = new BitSet();
		initialStates.set(states.nextSetBit(0));

		boolean changed = true;
		while (changed) {
			changed = false;
			actions.clear();
			for (int i = 0; i < model.getNumStates(); i++) {
				BitSet act = new BitSet();
				if (states.get(i)) {
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

	private List<BitSet> computeSCCs(NondetModel model) throws PrismException
	{
		SCCComputer sccc = SCCComputer.createSCCComputer(this, model);
		sccc.computeSCCs();
		return sccc.getSCCs();
	}

	private List<BitSet> translateStates(SubNondetModel model, List<BitSet> sccs)
	{
		List<BitSet> r = new ArrayList<BitSet>();
		for (int i = 0; i < sccs.size(); i++) {
			BitSet set = sccs.get(i);
			BitSet set2 = new BitSet();
			r.add(set2);
			for (int j = set.nextSetBit(0); j >= 0; j = set.nextSetBit(j + 1)) {
				set2.set(model.translateState(j));

			}
		}
		return r;
	}

	private boolean isMEC(BitSet b)
	{
		if (b.isEmpty())
			return false;

		int state = b.nextSetBit(0);
		while (state != -1) {
			boolean atLeastOneAction = false;
			for (int i = 0; i < model.getNumChoices(state); i++) {
				if (model.allSuccessorsInSet(state, i, b)) {
					atLeastOneAction = true;
				}
			}
			if (!atLeastOneAction) {
				return false;
			}
			state = b.nextSetBit(state + 1);
		}

		return true;
	}

	/**Described in paper of Maxi Weininger, SIMCEC
	 * Finds simple CECs in the MEC given to it, fixing the decisions of the minimizer according to L
	 */
	public List<BitSet> getSimBCECStates(BitSet mec, double[] L) {
		SMG smg = new SMG((SMG) model); //only call this method on SMGs, else it is not needed, since in MDPs all MECs are simple BCECs

		for (int s = 0; s<smg.getNumStates(); s++){
			if (smg.getPlayer(s)==2){
				for (int a = 0; a<smg.getNumChoices(s);a++){
					if (valueForStateActionPair(s,a,L,smg) > L[s]){
						smg.disableChoice(s,a);
					}
				}
			}
		}

		try {
			explicit.ECComputerDefault ec = (ECComputerDefault) ECComputer.createECComputer(this, smg);
			ec.computeMECStates(mec);
			return ec.getMECStates();
		} catch (PrismException e) {
			return null;
		}
	}

	public double getMinStayingValue(int s, BitSet mec, SMG smg, double[] vector){
		double result = 2;//too large, if 2 we know there is no staying L
		for (int a=0; a<smg.getNumChoices(s);a++){
			double Lsa = valueForStateActionPair(s,a,vector,smg);
			if(!leaves(s,a,mec) && Lsa<result){
				result = Lsa;
			}
		}
		return result;
	}

	public double getMaxLeavingValue(int s, BitSet mec, SMG smg, double[] vector){
		double result = 0;//no leaving action => value 0
		for (int a=0; a<smg.getNumChoices(s);a++){
			double Usa = valueForStateActionPair(s,a,vector,smg);
			if(!leaves(s,a,mec) && Usa>result){
				result = Usa;
			}
		}
		return result;
	}

	//Computes the U/L of a state action pair, given state, action, the vector U/L and the smg.
	public double valueForStateActionPair(int s, int a, double[] vect, SMG smg){
		Distribution distr = smg.trans.get(s).get(a);
		double result = 0.0;
		for (Map.Entry<Integer, Double> e : distr) {
			int s_prime = (Integer) e.getKey();
			double prob = (Double) e.getValue();
			result += prob * vect[s_prime];
		}
		return result;
	}

	public boolean leaves(int s, int a, BitSet mec){
		return !model.allSuccessorsInSet(s,a,mec);
	}
}
