//==============================================================================
//
//	Copyright (c) 2016-
//	Authors:
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

import common.FastUtils;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.IntArrayFIFOQueue;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntPriorityQueue;
import it.unimi.dsi.fastutil.objects.Object2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import jltl2ba.APElement;
import jltl2ba.APSet;
import jltl2dstar.NBA;
import jltl2dstar.NBA_State;
import prism.PrismLog;

import java.util.ArrayList;
import java.util.List;

/**
 * Construction and storage for the product
 * of the underlying labeled transition system
 * of a model (i.e., the edge relation defined
 * by getSuccessorsIterator()) and a non-deterministic
 * Büchi automaton (NBA).
 */
public class LTSNBAProduct extends Product<Model>
{
	/** A product state */
	private static class ProductState {
		/** The model state */
		private final int modelState;
		/** The automaton state */
		private final int automatonState;

		/** Constructor */
		public ProductState(int modelState, int automatonState)
		{
			this.modelState = modelState;
			this.automatonState = automatonState;
		}

		/** Get model state */
		public int getModelState()
		{
			return modelState;
		}

		/** Get automaton state */
		public int getAutomatonState()
		{
			return automatonState;
		}

		public String toString() {
			return "(" + getModelState() + "," + getAutomatonState() +")";
		}

		@Override
		public int hashCode()
		{
			final int prime = 31;
			int result = 1;
			result = prime * result + automatonState;
			result = prime * result + modelState;
			return result;
		}

		@Override
		public boolean equals(Object obj)
		{
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (!(obj instanceof ProductState))
				return false;
			ProductState other = (ProductState) obj;
			if (automatonState != other.automatonState)
				return false;
			if (modelState != other.modelState)
				return false;
			return true;
		}
	}

	/** Map state in product -> product state information */
	private final List<ProductState> productStates;
	/** Accepting Büchi states in the product */
	private final NatBitSet acceptingStates;

	/** Constructor for storing the product */
	private LTSNBAProduct(LTS productModel, Model originalModel, List<ProductState> productStates, NatBitSet acceptingStates)
	{
		super(productModel, originalModel);
		this.productStates = productStates;
		this.acceptingStates = acceptingStates;
	}

	@Override
	public int getModelState(int productState)
	{
		return productStates.get(productState).modelState;
	}

	@Override
	public int getAutomatonState(int productState)
	{
		return productStates.get(productState).automatonState;
	}

	/** Print the mapping between product state indizes and product states */
	public void printStateMapping(PrismLog log)
	{
		for (int i =0; i < productStates.size(); i++) {
			log.println(i + ": " + productStates.get(i) + (acceptingStates.contains(i) ? " !":""));
		}
	}

	/** Get the Büchi states for the product */
	public NatBitSet getAcceptingStates()
	{
		return acceptingStates;
	}

	/**
	 * Construct the product.
	 *
	 * @param model the model
	 * @param nba the nondeterministic Büchi automaton
	 * @param statesOfInterest the states in the model that serve as the starting point for the product
	 * @param labelBS vector of state sets for the atomic propositions L0, L1, ... in the automaton
	 */
	public static LTSNBAProduct doProduct(Model model, NBA nba, NatBitSet statesOfInterest, List<NatBitSet> labelBS)
	{
		// map state index in product automaton -> ProductState
		ArrayList<ProductState> productIdToProductState = new ArrayList<>();
		// map ProductState -> state index in product automaton
		Object2IntMap<ProductState> productStateToProductId = new Object2IntLinkedOpenHashMap<>();

		// storage for the product model
		LTSExplicit productModel = new LTSExplicit();
		// the accepting states in the product model
		NatBitSet acceptingStates = NatBitSets.set();

		// the stack of product ids that potentially have to be expanded
		IntPriorityQueue todo = new IntArrayFIFOQueue();

		// product state ids that have already been expanded
		NatBitSet expanded = NatBitSets.set();

		// Note: In contrast to the deterministic automaton product elsewhere in
		// PRISM, the automaton step delta(q,l(s)) is done when leaving the product
		// state, not when entering. As we are only interested in the infinite behaviour,
		// this doesn't change acceptance at all:
		//
		//     <s,q> -> <s',q'>  where s->s' in the model and q' in delta(q, l(s))
		//
		// Thus, the initial states in the product are <s, q_0> for s in the states of interest
		// of the model and q_0 the initial state of the NBA.

		NBA_State nbaStart = nba.getStartState();
		if (nbaStart == null) {
			// no start state = rejecting
			// to simplify treatment, add a dummy start state
			// with no outgoing transitions
			nbaStart = nba.newState();
			nba.setStartState(nbaStart);
		}

		// Note: As the NBA currently is guaranteed to have a single initial state,
		// the product has exactly one initial product state per model state of interest.
		// This allows the use of the "normal" projection back from the product result to
		// the original model. For multiple initial NBA states, it would be necessary
		// to aggregate the results over all initial NBA states when projecting back to
		// the original model.

		// for each model state of interest ...
		IntIterator iterator = FastUtils.iterator(statesOfInterest, model.getNumStates());
		while (iterator.hasNext()) {
			int modelState = iterator.nextInt();
			// ... construct a product state (modelState, nbaStart)
			ProductState p = new ProductState(modelState, nbaStart.getName());

			productIdToProductState.add(p);
			int id = productModel.addState();
			productModel.addInitialState(id);
			assert id == productIdToProductState.size() - 1;
			productStateToProductId.put(p, id);
			todo.enqueue(id);
			if (nbaStart.isFinal()) {
				acceptingStates.set(id);
			}
		}

		final APSet nbaAPSet = nba.getAPSet();

		while (!todo.isEmpty()) {
			int fromId = todo.dequeueInt();
			if (expanded.contains(fromId))
				continue;

			ProductState from = productIdToProductState.get(fromId);

			// construct edge label from the labeling of the model state
			APElement label = new APElement(nbaAPSet.size());
			for (int k = 0; k < nbaAPSet.size(); k++) {
				label.set(k, labelBS.get(Integer.parseInt(nbaAPSet.getAP(k).substring(1))).contains(from.getModelState()));
			}

			// the current state in the NBA
			NBA_State fromNBA = nba.get(from.getAutomatonState());
			// the successors in the NBA for the label corresponding to the model state
			NatBitSet nbaSuccessors = NatBitSets.asSet(fromNBA.getEdge(label));
			// for each successor of the model state ...
			IntIterator successorsIterator = model.getSuccessorsIterator(from.getModelState());
			while (successorsIterator.hasNext()) {
				int modelTo = successorsIterator.nextInt();
				// ... and NBA successor ...
				IntIterator nbaIterator = nbaSuccessors.iterator();
				while (nbaIterator.hasNext()) {
					int nbaSuccessor = nbaIterator.nextInt();
					// ... construct product state
					ProductState successor = new ProductState(modelTo, nbaSuccessor);
					int successorID = productStateToProductId.getOrDefault(successor, -1);
					if (successorID == -1) {
						// newly discovered, add as state to the product model
						productIdToProductState.add(successor);
						int newStateId = productModel.addState();
						assert newStateId == productIdToProductState.size() - 1;
						productStateToProductId.put(successor, newStateId);
						// mark as todo
						todo.enqueue(newStateId);
						boolean isAccepting = nba.get(nbaSuccessor).isFinal();
						if (isAccepting) {
							acceptingStates.set(newStateId);
						}
						successorID = newStateId;
					}

					// add the edge
					productModel.addEdge(fromId, successorID);
				}

				// fromId is fully expanded in the product
				expanded.set(fromId);
			}
		}

		return new LTSNBAProduct(productModel, model, productIdToProductState, acceptingStates);
	}
}
