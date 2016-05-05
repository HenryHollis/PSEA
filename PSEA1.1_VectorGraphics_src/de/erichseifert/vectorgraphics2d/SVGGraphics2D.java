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
package de.erichseifert.vectorgraphics2d;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Image;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RoundRectangle2D;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Map;

import javax.imageio.ImageIO;
import javax.xml.bind.DatatypeConverter;

import de.erichseifert.vectorgraphics2d.svg.SVGProcessor;
import de.erichseifert.vectorgraphics2d.util.DataUtils;
import de.erichseifert.vectorgraphics2d.util.GraphicsUtils;

/**
 * {@code Graphics2D} implementation that saves all operations to a string
 * in the <i>Scaled Vector Graphics</i> (SVG) format.
 */
@SuppressWarnings("restriction")
public class SVGGraphics2D extends ProcessingPipeline {
	private final Processor processor;

	/**
	 * Initializes a new VectorGraphics2D pipeline for translating Graphics2D
	 * commands to SVG data. The document dimensions must be specified as
	 * parameters.
	 * @param x Left offset.
	 * @param y Top offset
	 * @param width Width.
	 * @param height Height.
	 */
	public SVGGraphics2D(double x, double y, double width, double height) {
		super(x, y, width, height);
		processor = new SVGProcessor();
	}

	@Override
	protected Processor getProcessor() {
		return processor;
	}
}
