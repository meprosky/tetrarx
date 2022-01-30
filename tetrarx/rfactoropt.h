
extern DXRCL dxrcl;

inline float funcminrfactor(const nVector &point, void *tag)
{
	//Polyhedron probe(0.0f, LENEDGE, 0);

	//float f = dxrcl.rfactor_wcalcdxrcl();

	Polyhedron *p0;
	Polyhedron probe(0.0f, LENEDGE, 0);


	p0 = static_cast<Polyhedron*>(tag);
	//p0 = const_cast<Polyhedron*>(tag);

	int v1 = (*p0).fromvertex;
	int v2 = (*p0).tovertex;
	float rfactor0 = (*p0).rfactor;

	
	probe.connectcoaxandparallel(v1, *(*p0).topoly, v2);
	probe.rotate(probe.gc(v1), probe.gc(v1) - probe.gc(0), point.x[0]);
	probe.rotate(probe.gc(v1), probe.ye - probe.r0, point.x[1]);

	tetra[Ntetra] = probe;
	++Ntetra;
	//writetetra("temp.xyz");
	float rfactor = dxrcl.rfactor_wcalcdxrcl();
	--Ntetra;

	(*(*p0).topoly).a[v2].connected = 0;
    (*(*p0).topoly).a[v2].allowtoconnect = 1;

	return rfactor;

}


