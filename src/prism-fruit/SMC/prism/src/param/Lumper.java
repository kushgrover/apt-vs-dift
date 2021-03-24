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

import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.IntConsumer;

/**
 * Base class for signature-based bisimulation lumping.
 * Derived classes could implement specific types of bisimulation,
 * e.g. weak, strong, etc.
 *
 * @author Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
 * @see NullLumper
 * @see StrongLumper
 * @see WeakLumper
 */
abstract class Lumper {
	/**
	 * Kind of lumping to be used.
	 */
	enum BisimType {
		NULL,
		STRONG,
		WEAK
	}

	/** stores the original model to be lumped */
	protected MutablePMC origPmc;
	/** stores the lumped model */
	protected MutablePMC optPmc;
	/** stores the state partitioning */
	protected Partition partition;
	/** maintains a mapping from blocks to block numbers */
	protected Object2IntMap<IntSet> blockToNumber;
	/** list of all blocks */
	protected List<IntSet> blocks;
	/** maps states from original to containing states in lumped model */
	protected int[] originalToOptimised;

	/**
	 * Creates new lumper for given parametric Markov chain.
	 *
	 * @param origPmc Markov chain to build lumper for
	 */
	Lumper(MutablePMC origPmc)
	{
		this.origPmc = origPmc;
		partition = new Partition(origPmc);
		createInitialPartition();
		lump();
		mapBlocksToNumber();
		buildQuotient();
	}

	/**
	 * Creates a new lumper.
	 * This constructor needs to be present for technical reasons.
	 */
	protected Lumper()
	{

	}

	/**
	 * Stores combination of state reward and state time.
	 * This class is used subsume states with the same reward and time.
	 * If the time entry is irrelevant, it can be set to null.
	 *
	 * @author Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
	 */
	class RewardEntry
	{
		/** reward for this entry */
		final Function reward;
		/** time for this entry */
		final Function time;

		/**
		 * Creates a new reward entry.
		 *
		 * @param reward reward for this entry
		 * @param time time for this entry.
		 */
		RewardEntry(Function reward, Function time) {
			this.reward = reward;
			this.time = time;
		}

		@Override
		public int hashCode()
		{
			int hash = reward.hashCode() << 13;
			if (time != null) {
				hash |= time.hashCode();
			}
			return hash;
		}

		@Override
		public boolean equals(Object obj)
		{
			if (!(obj instanceof RewardEntry)) {
				return false;
			}

			RewardEntry other = (RewardEntry) obj;
			if (!this.reward.equals(other.reward)) {
				return false;
			}
			if ((this.time == null) != (other.time == null)) {
				return false;
			}
			if (this.time != null && !this.time.equals(other.time)) {
				return false;
			}
			return true;
		}
	}

	/**
	 * Creates an initial partitioning appropriate for analysis to perform.
	 * For reachability, takes care that only target-states / non-target
	 * are stored in one initial block, for reward analyses takes care
	 * that only states with same reward are subsumed, etc. For certain types
	 * of bisimulations, this is not enough, and the initial partition must be
	 * refined further before starting the main refinement process. This is
	 * the case for instance for weak bisimulation.
	 */
	protected void createInitialPartition()
	{
		IntSet oldBlock = partition.nextChangeableBlock();
		if (!origPmc.isUseRewards()) {
			IntSet safeStates = new IntOpenHashSet();
			IntSet unsafeStates = new IntOpenHashSet();
			oldBlock.forEach((IntConsumer) state -> {
				if (origPmc.isTargetState(state)) {
					unsafeStates.add(state);
				} else {
					safeStates.add(state);
				}
			});
			List<IntSet> newBlocks = new ArrayList<>();
			if (!safeStates.isEmpty()) {
				newBlocks.add(safeStates);
			}
			if (!unsafeStates.isEmpty()) {
				newBlocks.add(unsafeStates);
			}
			partition.addBlocks(newBlocks);
		} else {
			Map<RewardEntry, IntSet> rewardToStateMap = new HashMap<>();
			IntIterator iterator = oldBlock.iterator();
			while (iterator.hasNext()) {
				int state = iterator.nextInt();
				Function reward = origPmc.getReward(state);
				Function time = null;
				if (origPmc.isUseTime()) {
					time = origPmc.getTime(state);
				}
				RewardEntry entry = new RewardEntry(reward, time);
				IntSet block = rewardToStateMap.computeIfAbsent(entry, k -> new IntOpenHashSet());
				block.add(state);
			}
			List<IntSet> newBlocks = new ArrayList<>();
			for (IntSet block : rewardToStateMap.values()) {
				if (!block.isEmpty()) {
					newBlocks.add(block);
				}
			}
			partition.addBlocks(newBlocks);
		}
		partition.markAllBlocksAsNew();
	}

	/**
	 * Return the quotient of the original model.
	 * Before calling this method, {@code lump()} must have been called.
	 *
	 * @return quotient of original model
	 */
	MutablePMC getQuotient()
	{
		return optPmc;
	}

	/**
	 * Enumerates the blocks after refinement is finished.
	 * This is done as a preparation to create the quotient and the
	 * mapping from original to quotient states.
	 */
	protected void mapBlocksToNumber()
	{
		blocks = partition.getAllBlocks();
		blockToNumber = new Object2IntLinkedOpenHashMap<>();
		originalToOptimised = new int[origPmc.getNumStates()];
		int nextBlockNumber = 0;
		for (IntSet block : blocks) {
			blockToNumber.put(block, nextBlockNumber);
			IntIterator iterator = block.iterator();
			while (iterator.hasNext()) {
				originalToOptimised[iterator.nextInt()] = nextBlockNumber;
			}
			nextBlockNumber++;
		}
	}

	/**
	 * Returns a mapping from original model states to quotient states.
	 * Thus, {@code result[state]} is the state of the quotient to which
	 * {@code state} belongs.
	 *
	 * @return mapping from original model states to quotient states
	 */
	int[] getOriginalToOptimised()
	{
		return originalToOptimised;
	}

	/**
	 * Perform the lumping process.
	 * Refines blocks from the list of blocks to refine, it is clear that
	 * further tries to refine blocks will not lead to the creation of any
	 * new blocks.
	 */
	protected void lump()
	{
		ArrayList<IntSet> newBlocks = new ArrayList<>();
		while (partition.mayChange()) {
			newBlocks.clear();
			IntSet oldBlock = partition.nextChangeableBlock();
			refineBlock(oldBlock, newBlocks);
			partition.addBlocks(newBlocks);
		}
	}

	/**
	 * Refine a given block.
	 * The exact way of how this is done depends on the type of bisimulation used.
	 *
	 * @param oldBlock block to refine
	 * @param newBlocks list of new blocks generated
	 */
	protected abstract void refineBlock(IntSet oldBlock, List<IntSet> newBlocks);

	/**
	 * Build the quotients after refinement has finished.
	 * Each state of the quotient will be a block of states from the original
	 * model. The exact way of how transitions between states are computed
	 * depends on the type of bisimulation used.
	 */
	protected abstract void buildQuotient();
}
