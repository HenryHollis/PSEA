/*
 * VectorGraphics2D: Vector export for Java(R) Graphics2D
 *
 * (C) Copyright 2010-2015 Erich Seifert <dev[at]erichseifert.de>
 *
 * This file is part of VectorGraphics2D.
 *
 * VectorGraphics2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * VectorGraphics2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with VectorGraphics2D.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.erichseifert.vectorgraphics2d.pdf;

import java.util.LinkedHashMap;
import java.util.Map;

public class PDFObject {
	public final int id;
	public final int version;
	public final Map<String, Object> dict;
	public final Payload payload;

	public PDFObject(int id, int version, Map<String, Object> dict, Payload payload) {
		this.dict = new LinkedHashMap<String, Object>();
		this.id = id;
		this.version = version;
		this.payload = payload;
		if (dict != null) {
			this.dict.putAll(dict);
		}
	}
}

