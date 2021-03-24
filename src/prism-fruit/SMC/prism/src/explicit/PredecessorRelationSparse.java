//==============================================================================
//
//	Copyright (c) 2014-
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

import de.tum.in.naturals.set.NatBitSet;
import it.unimi.dsi.fastutil.ints.Int2ObjectFunction;
import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntLists;

import java.util.function.IntConsumer;

public class PredecessorRelationSparse implements PredecessorRelation
{
	private final Model model;
	Int2ObjectFunction<IntCollection> pre;

	public PredecessorRelationSparse(Model model, NatBitSet subset)
	{
		this.model = model;
		pre = new Int2ObjectLinkedOpenHashMap<>(subset.size());
		subset.forEach((IntConsumer) state -> {
			IntIterator successorIterator = model.getSuccessorsIterator(state);
			successorIterator.forEachRemaining((IntConsumer) successor -> {
				IntCollection predecessors = pre.get(successor);
				if (predecessors == null) {
					predecessors = new IntArrayList();
					pre.put(successor, predecessors);
				}
				assert !predecessors.contains(state);
				predecessors.add(state);
			});
		});
		pre.defaultReturnValue(IntLists.EMPTY_LIST);
	}

	@Override public IntIterable getPre(int s)
	{
		return pre.get(s);
	}

	@Override public IntIterator getPredecessorsIterator(int s)
	{
		return pre.get(s).iterator();
	}

	@Override public Model getModel()
	{
		return model;
	}
}
