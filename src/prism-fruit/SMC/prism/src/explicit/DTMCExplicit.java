//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import prism.PrismLog;

import java.util.Iterator;

/**
 * Base class for explicit-state representations of a DTMC.
 */
public abstract class DTMCExplicit extends ModelExplicit implements DTMC
{
	// Accessors (for Model)

	@Override
	protected void exportTransitionsToDotFile(int i, PrismLog out)
	{
		Iterator<Int2DoubleMap.Entry> iter = getTransitionsIterator(i);
		while (iter.hasNext()) {
			Int2DoubleMap.Entry e = iter.next();
			out.print(i + " -> " + e.getIntKey() + " [ label=\"");
			out.print(e.getDoubleValue() + "\" ];\n");
		}
	}

	// Accessors (for DTMC)

	@Override
	public String toString()
	{
		boolean first = true;
		StringBuilder s = new StringBuilder("trans: [ ");
		for (int i = 0; i < numStates; i++) {
			if (first) {
				first = false;
			} else {
				s.append(", ");
			}
			s.append(i).append(": {");
			if (getNumTransitions(i) != 0) {
				Iterator<Int2DoubleMap.Entry> transitionsIterator = getTransitionsIterator(i);
				while (true) {
					Int2DoubleMap.Entry trans = transitionsIterator.next();
					s.append(trans.getIntKey()).append("=").append(trans.getDoubleValue());
					if (transitionsIterator.hasNext()) {
						s.append(", ");
					} else {
						break;
					}
				}
			}
			s.append("}");
		}
		s.append(" ]");
		return s.toString();
	}
}
