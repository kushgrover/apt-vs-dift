package org.ojalgo.matrix.store;

import org.ojalgo.array.SparseArray;

public class SparseStore2
{
	public static SparseStore<Double> makePrimitive(int rowsCount, int columnsCount) {
		return new SparseStore<>(PrimitiveDenseStore.FACTORY, rowsCount, columnsCount, SparseArray.makePrimitive((long) rowsCount * columnsCount));
	}
}
