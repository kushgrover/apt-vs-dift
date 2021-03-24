//==============================================================================
//
//	Copyright (c) 2013-
//	Authors:
//	* Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
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

package param;

import it.unimi.dsi.fastutil.ints.Int2ObjectFunction;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayFIFOQueue;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntListIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntPriorityQueue;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * Weak bisimulation lumper.
 * Notice that weak bisimulation is only valid for unbounded reachability,
 * but must not be used for expected accumulated rewards or long-run
 * average rewards.
 *
 * @author Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
 */
final class WeakLumper extends Lumper {

	/**
	 * Construct a new weak bisimulation lumper.
	 *
	 * @param origPmc Markov chain to construct lumper for
	 */
	WeakLumper(MutablePMC origPmc) {
		super(origPmc);
	}

	/**
	 * Construct the weak bisimulation signature of given block.
	 * The signature is a mapping of blocks to the probability to move
	 * from the given state to any state of the block under the condition
	 * that the block is left. Thus, it is only defined for states which have
	 * a non-zero probability of leaving the block in one step. The function
	 * returns {@code null} for states ("silent" states) for which this is
	 * not the case.
	 *
	 * @param state state to compute signature of
	 * @return signature of this state
	 */
	private Map<IntSet, Function> stateSignature(int state, IntSet ownClass)
	{
		Map<IntSet, Function> signature = new HashMap<>();
		IntListIterator toStateIter = origPmc.transitionTargets.get(state).listIterator();
		ListIterator<Function> toProbIter = origPmc.transitionProbs.get(state).listIterator();

		/* compute probability to remain in block in one step */
		Function slProb = origPmc.getFunctionFactory().getZero();
		while (toStateIter.hasNext()) {
			int toState = toStateIter.nextInt();
			Function toStateProb = toProbIter.next();
			if (ownClass.contains(toState)) {
				slProb = slProb.add(toStateProb);
			}
		}
		/* for states which cannot leave their block directly, return {@code null} */
		if (slProb.equals(origPmc.getFunctionFactory().getOne())) {
			return null;
		}
		/* 1 / (1 - slProb) */
		Function star = slProb.star();

		toStateIter = origPmc.transitionTargets.get(state).listIterator();
		toProbIter = origPmc.transitionProbs.get(state).listIterator();
		while (toStateIter.hasNext()) {
			int toState = toStateIter.nextInt();
			Function toStateProb = toProbIter.next();
			IntSet toBlock = partition.getStateBlock(toState);
			if (ownClass != toBlock) {
				toStateProb = star.multiply(toStateProb);
				Function toBlockProb = signature.get(toBlock);
				if (toBlockProb == null) {
					toBlockProb = origPmc.getFunctionFactory().getZero();
				}
				toBlockProb = toBlockProb.add(toStateProb);
				signature.put(toBlock, toBlockProb);
			}
		}
		return signature;
	}

	/**
	 * Refines a given block to a list of new blocks for weak bisimulation.
	 * New blocks are as follows: some of the new blocks consist of the
	 * states which can leave their block ("non-silent" states) and which
	 * have the same signature. In addition, such a block contains the states
	 * which cannot leave their block in one step ("silent" states) and which
	 * can only reach states of this particular new block. Other blocks
	 * consist of silent states which can reach more than one particular
	 * signature block. For these kind of blocks, we have a new block for
	 * each combination of new blocks they might reach. For instance, if
	 * there are new blocks (based on a signature) A,B,C, we add blocks
	 * {A,B},{B,C} and {A,B,C}, containing silent states which can reach
	 * A and B / B and C / A, B and C.
	 *
	 * @param oldBlock block to refine
	 * @param newBlocks list of new blocks generated
	 */
	@Override
	protected void refineBlock(IntSet oldBlock, List<IntSet> newBlocks) {
		IntList nonSilent = new IntArrayList(oldBlock.size());
		IntSet silent = new IntOpenHashSet();
		Map<Map<IntSet, Function>, IntSet> signatures = new HashMap<>();
		Int2ObjectFunction<IntSet> stateToBlock = new Int2ObjectOpenHashMap<>();
		/* compute signatures of states of old block and divide into silent/
		 * nonsilent states. Silent states are states which cannot leave the
		 * block in one step. */
		oldBlock.forEach((int state) -> {
			Map<IntSet, Function> signature = stateSignature(state, oldBlock);
			if (signature != null) {
				nonSilent.add(state);
				IntSet newBlock = signatures.computeIfAbsent(signature, k -> new IntOpenHashSet());
				newBlock.add(state);
				stateToBlock.put(state, newBlock);
			} else {
				silent.add(state);
			}
		});

		/* non-silent states reach only the new block they are contained in */
		Int2ObjectMap<Set<IntSet>> reachWhichBlocks = new Int2ObjectOpenHashMap<>();
		oldBlock.forEach((int state) -> {
			Set<IntSet> predReachBlocks = new HashSet<>();
			if (!silent.contains(state)) {
				predReachBlocks.add(stateToBlock.get(state));
			}
			reachWhichBlocks.put(state, predReachBlocks);
		});

		/* collect all silent states which can reach a particular
		 * non-silent state by performing a backwards depth-first search.
		 * Mark silent states one comes across with the block of the
		 * non-silent state. We can already stop the search if we know
		 * that the state has previously been visited from another
		 * state from the same block. */
		nonSilent.forEach((int state) -> {
			IntSet block = stateToBlock.get(state);
			IntPriorityQueue stack = new IntArrayFIFOQueue();
			stack.enqueue(state);
			while (!stack.isEmpty()) {
				int stackState = stack.dequeueInt();
				origPmc.incoming.get(stackState).forEach((int predState) -> {
					Set<IntSet> predReachBlocks = reachWhichBlocks.get(predState);
					if (oldBlock.contains(predState) && silent.contains(predState) && !predReachBlocks.contains(block)) {
						predReachBlocks.add(block);
						stack.enqueue(predState);
					}
				});
			}
		});

		/* compute new blocks, add the nonempty ones to list of new blocks */
		Map<Set<IntSet>, IntSet> remap = new HashMap<>();
		for (Int2ObjectMap.Entry<Set<IntSet>> entry : reachWhichBlocks.int2ObjectEntrySet()) {
			IntSet sigStates = remap.computeIfAbsent(entry.getValue(), k -> new IntOpenHashSet());
			sigStates.add(entry.getIntKey());
		}
		for (IntSet block : remap.values()) {
			if (!block.isEmpty()) {
				newBlocks.add(block);
			}
		}
	}

	/**
	 * Build the weak bisimulation quotient from the blocks computed.
	 * Transition probabilities are basically based on the weak bisimulation
	 * signature. However, we must take care that we use a non-silent state
	 * to compute transition probabilties from. Also, states which can never
	 * leave their block after an arbitrary number of steps ("divergent"
	 * states) must lead to adding a self loop in their containing block.
	 */
	@Override
	protected void buildQuotient() {
		optPmc = new MutablePMC(origPmc.getFunctionFactory(), blocks.size(), origPmc.isUseRewards(), false);
		for (int newState = 0; newState < blocks.size(); newState++) {
			Map<IntSet, Function> signature = null;
			int oldState = -1;
			IntSet fromBlock = blocks.get(newState);
			IntIterator iter = fromBlock.iterator();
			while (iter.hasNext()) {
				oldState = iter.nextInt();
				signature = stateSignature(oldState, fromBlock);
				if (signature != null) {
					break;
				}
			}
			if (signature == null) {
				optPmc.addTransition(newState, newState, origPmc.getFunctionFactory().getOne());
			} else {
				for (Entry<IntSet, Function> entry : signature.entrySet()) {
					optPmc.addTransition(newState, blockToNumber.getInt(entry.getKey()), entry.getValue());
				}
			}
			if (origPmc.isUseRewards()) {
				optPmc.setReward(newState, origPmc.getReward(oldState));
			} else {
				optPmc.setTargetState(newState, origPmc.isTargetState(oldState));
			}
		}
	}

	/**
	 * Creates an initial partitioning.
	 * This function is based on the of the {@code Lumper} class. However,
	 * for the weak bisimulation lumping to work correctly, for each block
	 * of the initial partitioning, we have to split off "divergent" states.
	 * Divergent states are states which can never leave their block, after
	 * any number of steps.
	 */
	@Override
	protected void createInitialPartition() {
		super.createInitialPartition();
		List<IntSet> newBlocks = new ArrayList<>();
		while (partition.mayChange()) {
			IntSet oldBlock = partition.nextChangeableBlock();
			IntSet leaveSet = new IntOpenHashSet();
			IntArrayList directLeaving = new IntArrayList();
			oldBlock.forEach((int state) -> {
				IntIterator iterator = origPmc.transitionTargets.get(state).iterator();
				while (iterator.hasNext()) {
					int toState = iterator.nextInt();
					if (!oldBlock.contains(toState)) {
						leaveSet.add(state);
						directLeaving.add(state);
						break;
					}
				}
			});
			directLeaving.forEach((int state) -> {
				IntArrayFIFOQueue stack = new IntArrayFIFOQueue();
				stack.enqueue(state);
				while (!stack.isEmpty()) {
					int leaving = stack.dequeueInt();
					origPmc.incoming.get(leaving).forEach((int inState) -> {
						if (oldBlock.contains(inState) && !leaveSet.contains(inState)) {
							leaveSet.add(inState);
							stack.enqueue(inState);
						}
					});
				}
			});
			if (!leaveSet.isEmpty()) {
				newBlocks.add(leaveSet);
			}
			IntSet staying = new IntOpenHashSet(oldBlock);
			staying.removeAll(leaveSet);
			if (!staying.isEmpty()) {
				newBlocks.add(staying);
			}
		}
		partition.addBlocks(newBlocks);
		partition.markAllBlocksAsNew();
	}
}
