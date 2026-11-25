import sys, importlib.util, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"../demo_2d"); sys.path.insert(0,"../src")
import cases3d as K3, surface3d as S3, gllbase
# import run_demo WITHOUT running main()
spec=importlib.util.spec_from_file_location("rd","run_demo.py")
rd=importlib.util.module_from_spec(spec)
sys.modules["rd"]=rd
src=open("run_demo.py").read().replace('if __name__ == "__main__":\n    main()','')
exec(compile(src,"run_demo.py","exec"),rd.__dict__)

def area(tris):
    return sum(0.5*np.linalg.norm(np.cross(np.array(b)-np.array(a),np.array(c)-np.array(a)))
               for (a,b,c) in tris)
print("%-10s %11s %11s %10s %8s"%("case","extracted","reference","deficit","%"))
for key in ("tgv","turbulent"):
    c=K3.CASES[key]; N=K3.case_order(key)
    xn,_=gllbase.gll(N); R,S,T=np.meshgrid(xn,xn,xn,indexing="ij")
    F=c["f"](R,S,T)-c["level"]
    A=area(S3.extract_surface(N,F,0.0).tris)
    Aref=area(rd.analytic_tris(key,ngrid=100))
    print("%-10s %11.3f %11.3f %10.3f %7.1f%%"%(key,A,Aref,Aref-A,100*(Aref-A)/Aref))
