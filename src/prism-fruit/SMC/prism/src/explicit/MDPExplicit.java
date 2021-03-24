//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Christian von Essen <christian.vonessen@imag.fr> (Verimag, Grenoble)
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

/**
 * Base class for explicit-state representations of an MDP.
 */
public abstract class MDPExplicit extends ModelExplicit implements MDP
{
	// Accessors (for Model)

	@Override
	public void exportTransitionsToDotFile(int i, PrismLog out)
	{
		int j, numChoices;
		String nij;
		Object action;
		numChoices = getNumChoices(i);
		for (j = 0; j < numChoices; j++) {
			action = getAction(i, j);
			nij = "n" + i + "_" + j;
			out.print(i + " -> " + nij + " [ arrowhead=none,label=\"" + j);
			if (action != null)
				out.print(":" + action);
			out.print("\" ];\n");
			out.print(nij + " [ shape=point,width=0.1,height=0.1,label=\"\" ];\n");
			for (Int2DoubleMap.Entry e : getTransitions(i, j)) {
				out.print(nij + " -> " + e.getIntKey() + " [ label=\"" + e.getDoubleValue() + "\" ];\n");
			}
		}
	}
}
