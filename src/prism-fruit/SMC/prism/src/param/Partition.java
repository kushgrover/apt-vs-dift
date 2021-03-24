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

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.function.IntConsumer;

/**
 * Stores a given partitioning of the states of a Markov model.
 * This class is intended to be used in combination with signature-based
 * partitioning refinement. In addition to storing the states, it can also
 * compute which blocks become subject to change due to changes in other
 * blocks.
 *
 * @author Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
 */
final class Partition {
	/**
	 * Comparator class comparing integer hash sets according to their size.
	 *
	 * @author Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
	 */
	final class HashSetSizeComparator implements Comparator<IntSet> {
		@Override
		public int compare(IntSet o1, IntSet o2) {
			int size1 = o1.size();
			int size2 = o2.size();
			// TODO should actually be other way round, but slower then... ??
			if (size1 < size2) {
				return 1;
			} else if (size1 == size2) {
				return 0;
			} else {
				return -1;
			}
		}

		@Override
		public boolean equals(Object obj) {
			return obj instanceof HashSetSizeComparator;
		}

		@Override
		public int hashCode() {
			return 0;
		}
	}

	/** parametric Markov chain this partition is for */
	private MutablePMC pmc;
	/** all blocks of this partitioning */
	private Set<IntSet> blocks;
	/** maps states to the block they are contained in */
	private List<IntSet> stateToBlock;
	/** list of blocks which might need to be refined */
	private PriorityQueue<IntSet> mayChange;
	/** hash set of blocks which might need to be refined */
	private Set<IntSet> mayChangeHash;
	/** next block to refine */
	private IntSet nextBlock;

	/**
	 * Creates new partitioning for given parametric Markov chain.
	 * Initially, all states will be in the same block.
	 *
	 * @param pmc parametric Markov chain to create partitioning of
	 */
	Partition(MutablePMC pmc)
	{
		this.pmc = pmc;
		blocks = new HashSet<IntSet>();
		stateToBlock = new ArrayList<IntSet>(pmc.getNumStates());
		IntSet initialBlock = new IntOpenHashSet(1);
		mayChange = new PriorityQueue<IntSet>(11, new HashSetSizeComparator());
		mayChangeHash = new HashSet<IntSet>();
		for (int state = 0; state < pmc.getNumStates(); state++) {
			initialBlock.add(state);
			stateToBlock.add(initialBlock);
		}
		blocks.add(initialBlock);
		mayChange.add(initialBlock);
		mayChangeHash.add(initialBlock);
	}

	/**
	 * Obtain the next changeable block, according to their priority (size).
	 *
	 * @return next changeable block
	 */
	IntSet nextChangeableBlock()
	{
		nextBlock = mayChange.poll();
		mayChangeHash.remove(nextBlock);
		blocks.remove(nextBlock);
		return nextBlock;
	}

	/**
	 * Add a list of new blocks to list of blocks.
	 * Also marks existing blocks as subject to change in case they have
	 * transitions into these new blocks, and thus might have to be split.
	 *
	 * @param newBlocks list of blocks to add
	 */
	void addBlocks(List<IntSet> newBlocks)
	{
		blocks.addAll(newBlocks);

		for (IntSet block : newBlocks) {
			block.forEach((IntConsumer) state -> stateToBlock.set(state, block));
			if (block.equals(nextBlock)) {
				return;
			}
			mayChange.add(block);
			mayChangeHash.add(block);
		}
		for (IntSet block : newBlocks) {
			block.forEach((IntConsumer) state -> {
				pmc.incoming.get(state).forEach((IntConsumer) predec -> {
					IntSet predecBlock = stateToBlock.get(predec);
					if (!nextBlock.contains(predec) && !mayChangeHash.contains(predecBlock) && (predecBlock.size() > 1)) {
						mayChange.add(predecBlock);
						mayChangeHash.add(predecBlock);
					}
				});
			});
		}
	}

	/**
	 * Checks whether there are blocks remaining which might need refinement.
	 *
	 * @return true iff blocks are remaining which might need refinement
	 */
	boolean mayChange()
	{
		return !mayChange.isEmpty();
	}

	/**
	 * Obtain a list of all blocks of the partition.
	 *
	 * @return list of all blocks of the partition
	 */
	ArrayList<IntSet> getAllBlocks()
	{
		ArrayList<IntSet> allBlocks = new ArrayList<IntSet>(blocks.size());
		for (IntSet block : blocks) {
			allBlocks.add(block);
		}

		return allBlocks;
	}

	/**
	 * Get the block in which a given state is contained.
	 *
	 * @param state state to search containing block of
	 * @return block state is contained in
	 */
	IntSet getStateBlock(int state)
	{
		return stateToBlock.get(state);
	}

	/**
	 * Marks all blocks as being new.
	 */
	void markAllBlocksAsNew()
	{
		mayChange.clear();
		mayChangeHash.clear();
		for (IntSet block : blocks) {
			mayChange.add(block);
			mayChangeHash.add(block);
		}
	}
}
