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
import it.unimi.dsi.fastutil.ints.Int2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntStack;

import java.util.ArrayList;
import java.util.List;

public final class Tarjan
{
	private final IntStack stack;
	private final NatBitSet onStack;
	private final TarjanModel model;
	private final SCCConsumer sccConsumer;
	private final NatBitSet todo;
	private final Int2IntMap lowLinks;
	private int index = 0;
	private boolean doneEarly = false;

	private Tarjan(TarjanModel model, SCCConsumer sccConsumer, NatBitSet states)
	{
		this.model = model;
		this.sccConsumer = sccConsumer;
		this.todo = states;
		stack = new IntArrayList();
		onStack = NatBitSets.set();
		lowLinks = new Int2IntLinkedOpenHashMap();
		lowLinks.defaultReturnValue(-1);
	}

	public static boolean isBSCC(NatBitSet scc, TarjanModel model)
	{
		IntIterator sccStates = scc.iterator();
		while (sccStates.hasNext()) {
			IntIterator successorsIterator = model.successors(sccStates.nextInt());
			while (successorsIterator.hasNext()) {
				if (!scc.contains(successorsIterator.nextInt())) {
					return false;
				}
			}
		}
		return true;
	}

	public static void runTarjan(TarjanModel model, SCCConsumer sccConsumer, NatBitSet states)
	{
		getTarjan(model, sccConsumer, states.clone()).run();
	}

	public static Tarjan getTarjan(TarjanModel model, SCCConsumer sccConsumer, NatBitSet states)
	{
		return new Tarjan(model, sccConsumer, states.clone());
	}

	public static ECComputerFast.DecompositionResult getSCCs(TarjanModel model, NatBitSet states)
	{
		List<NatBitSet> sccs = new ArrayList<>();
		List<NatBitSet> bsccs = new ArrayList<>();
		SCCConsumer consumer = scc -> {
			if (isBSCC(scc, model)) {
				bsccs.add(scc);
			} else {
				sccs.add(scc);
			}
			return true;
		};
		runTarjan(model, consumer, states);
		return new ECComputerFast.DecompositionResult(sccs, bsccs);
	}

	public static List<NatBitSet> getBSCCs(TarjanModel model, NatBitSet states)
	{
		List<NatBitSet> bsccs = new ArrayList<>();
		SCCConsumer consumer = scc -> {
			if (isBSCC(scc, model)) {
				bsccs.add(scc);
			}
			return true;
		};
		runTarjan(model, consumer, states);
		return bsccs;
	}

	boolean isTodo(int state)
	{
		return todo.contains(state);
	}

	boolean isOnStack(int state)
	{
		return onStack.contains(state);
	}

	public void run()
	{
		while (!todo.isEmpty()) {
			runFrom(todo.firstInt());
		}
	}

	// TODO: Replace recursive calls by while-loop and an int call stack.
	private void runFrom(int state)
	{
		assert !onStack.contains(state);

		int stateIndex = index;
		lowLinks.put(state, stateIndex);
		stack.push(state);
		onStack.set(state);
		todo.clear(state);
		index += 1;

		boolean updatedLowLink = false;

		// foreach loop has boxing
		IntIterator iterator = model.successors(state);
		while (iterator.hasNext()) {
			int successor = iterator.nextInt();
			if (successor == state) {
				// No need to process self-loops
				continue;
			}
			if (todo.contains(successor)) {
				// Successor was not processed, recurse
				runFrom(successor);

				int successorLowLink = lowLinks.get(successor);
				if (successorLowLink == -1) {
					// This successor is the root a fully discovered SCC or a deadlock state
					if (doneEarly) {
						return;
					}
					continue;
				}
				assert lowLinks.containsKey(state);
				if (successorLowLink < lowLinks.get(state)) {
					// Set the link of the node to the min of it and the successor's link
					lowLinks.put(state, successorLowLink);
					updatedLowLink = true;
				}
			} else if (onStack.contains(successor)) {
				// Successor is still on the stack, i.e. not fully explored and we found a link to it, hence the link of this state is <= than the
				// successors link
				assert lowLinks.containsKey(successor);
				int successorLowLink = lowLinks.get(successor);

				assert lowLinks.containsKey(state);
				if (successorLowLink < lowLinks.get(state)) {
					lowLinks.put(state, successorLowLink);
					updatedLowLink = true;
				}
			} else {
				// Successor was fully explored
				assert !lowLinks.containsKey(successor);
			}
		}

		assert !updatedLowLink == (lowLinks.get(state) == stateIndex);
		if (!updatedLowLink) {
			assert !stack.isEmpty() && onStack.contains(state);
			int sccState = popState();
			NatBitSet scc;
			if (sccState == state) {
				scc = NatBitSets.singleton(sccState);
			} else {
				scc = NatBitSets.set();
				scc.set(sccState);
				while (sccState != state) {
					sccState = popState();
					scc.set(sccState);
				}
			}

			if (!sccConsumer.accept(scc)) {
				doneEarly = true;
			}
		}
	}

	private int popState() {
		int state = stack.popInt();
		assert onStack.contains(state) && lowLinks.containsKey(state);
		onStack.clear(state);
		lowLinks.remove(state);
		return state;
	}

	@FunctionalInterface
	public interface SCCConsumer
	{
		/**
		 * Gets informed about each found SCC. If false is returned, the SCC search stops. Is called in order of detection, i.e. in reverse topological order.
		 */
		boolean accept(NatBitSet scc);
	}

	@FunctionalInterface
	public interface TarjanModel
	{
		IntIterator successors(int state);
	}

	public static final class NondetMecTarjanModel implements TarjanModel
	{
		private final NondetModel model;
		private final MEC mec;

		public NondetMecTarjanModel(NondetModel model, MEC mec)
		{
			this.model = model;
			this.mec = mec;
		}

		@Override public IntIterator successors(int state)
		{
			NatBitSet stateActions = mec.actions.get(state);
			if (stateActions.size() == 1) {
				return model.getSuccessors(state, stateActions.firstInt()).iterator();
			}

			IntIterator[] successors = new IntIterator[stateActions.size()];
			int index = 0;
			IntIterator actionIterator = stateActions.iterator();
			while (actionIterator.hasNext()) {
				successors[index] = model.getSuccessors(state, actionIterator.nextInt()).iterator();
				index += 1;
			}
			return IntIterators.concat(successors);
		}
	}

	public static final class DefaultTarjanModel implements TarjanModel
	{
		private final Model model;

		public DefaultTarjanModel(Model model)
		{
			this.model = model;
		}

		@Override public IntIterator successors(int state)
		{
			return model.getSuccessorsIterator(state);
		}
	}
}
