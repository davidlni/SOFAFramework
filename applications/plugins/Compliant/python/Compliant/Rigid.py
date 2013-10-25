import Sofa

from subprocess import Popen, PIPE


# small helper
def concat(x):
        return ' '.join(map(str, x))


class Frame:
        def __init__(self):
                self.translation = [0, 0, 0]
                self.rotation = [0, 0, 0, 1]
                
        def mstate(self, parent, **args):
                return parent.createObject('MechanicalObject', 
                                           template = 'Rigid',
                                           translation = concat(self.translation),
                                           rotation = concat(self.rotation),
                                           **args)
                
        def __str__(self):
                return concat(self.translation) + ' ' + concat(self.rotation) 
               
        def read(self, str):
                num = map(float, str.split())
                self.translation = num[:3]
                self.rotation = num[3:]
                return self

def mesh_offset( mesh, path = "" ):
        str = subprocess.Popen("GenerateRigid " + mesh ).stdout.read()
        

class MassInfo:
        pass

# density is kg/m^3
def generate_rigid(filename, density = 1000.0):
        cmd = Sofa.build_dir() + '/bin/GenerateRigid'
        args = filename
        output = Popen([cmd, args], stdout=PIPE)
        line = output.stdout.read().split('\n')
        start = 2
        
        # print line 
        
        mass = float( line[start].split(' ')[1] )
        volm = float( line[start + 1].split(' ')[1] )
        inrt = map(float, line[start + 2].split(' ')[1:] )
        com = map(float, line[start + 2].split(' ')[1:] )
        
        # TODO extract principal axes basis if needed
        # or at least say that we scred up
        
        res = MassInfo()

        # by default, GenerateRigid assumes 1000 kg/m^3 already
        res.mass = (density / 1000.0) * mass
        
        res.inertia = [inrt[0], inrt[3 + 1], inrt[6 + 2] ]
        res.inertia = [(density / 1000.0) * x for x in res.inertia]
        res.com = com

        return res

# TODO provide synchronized members with sofa states
class Body:
        
        def __init__(self, name = "unnamed"):
                self.name = name
                self.collision = None # collision mesh
                self.visual = None    # visual mesh
                self.dofs = Frame()
                self.mass = 1
                self.inertia = [1, 1, 1] 
                self.template = 'Rigid'
                self.color = [1, 1, 1]
                # TODO more if needed (scale, color)
                
        def mass_from_mesh(self, name, density = 1000.0):
                info = generate_rigid(name, density)
                self.mass = info.mass
                self.inertia = info.inertia
                
                # TODO handle com/principal basis

        def insert(self, node):
                res = node.createChild( self.name )

                dofs = self.dofs.mstate(res, name = 'dofs' )
                
                mass = res.createObject('RigidMass', 
                                        template = self.template, 
                                        name = 'mass', 
                                        mass = self.mass, 
                                        inertia = concat(self.inertia))

                
                # visual
                if self.visual != None:
                        visual_template = 'ExtVec3f'
                        
                        visual = res.createChild( 'visual' )
                        ogl = visual.createObject('OglModel', 
                                                  template = visual_template, 
                                                  name='mesh', 
                                                  fileMesh = self.visual, 
                                                  color = concat(self.color), 
                                                  scale3d='1 1 1')
                        
                        visual_map = visual.createObject('RigidMapping', 
                                                         template = self.template + ', ' + visual_template, 
                                                         input = '@../')
                
                if self.collision != None:
                        # collision
                        collision = res.createChild('collision')
                
                        # TODO lol
                
                return res


class Joint:

        def __init__(self, name = 'joint'):
                self.dofs = [0, 0, 0, 0, 0, 0]
                self.body = []
                self.offset = []
                self.damping = 0
                self.name = name
                self.compliance = 0
                
        def append(self, node, offset = None):
                self.body.append(node)
                self.offset.append(offset)
                self.name = self.name + '-' + node.name
        
        class Node:
                pass
        
        def insert(self, parent, add_compliance = True):
                # build input data for multimapping
                input = []
                for b, o in zip(self.body, self.offset):
                        if o is None:
                                input.append( '@' + b.name + '/dofs' )
                        else:
                                joint = b.createChild( self.name + '-offset' )
                                
                                joint.createObject('MechanicalObject', 
                                                   template = 'Rigid', 
                                                   name = 'dofs' )
                                
                                joint.createObject('AssembledRigidRigidMapping', 
                                                   template = "Rigid,Rigid",
                                                   source = '0 ' + str( o ) )
                                
                                input.append( '@' + b.name + '/' + joint.name + '/dofs' )
                                
                # now for the joint dofs
                node = parent.createChild(self.name)
                
                dofs = node.createObject('MechanicalObject', 
                                         template = 'Vec6d', 
                                         name = 'dofs', 
                                         position = '0 0 0 0 0 0' )
                
                # TODO handle damping
                mask = [ (1 - d) for d in self.dofs ]
                
                map = node.createObject('RigidJointMultiMapping',
                                        name = 'mapping', 
                                        template = 'Rigid,Vec6d', 
                                        input = concat(input),
                                        output = '@dofs',
                                        dofs = concat(mask),
                                        pairs = "0 0")
                
                if add_compliance:
                        compliance = node.createObject('UniformCompliance',
                                                       name = 'compliance',
                                                       template = 'Vec6d',
                                                       compliance = self.compliance)
                        stab = node.createObject('Stabilization')

                # for some reason return node is unable to lookup for
                # children using getChild() so in the meantime...
                res = Joint.Node()
                
                res.node = node
                # res.compliance = compliance
                
                return res

class SphericalJoint(Joint):

        def __init__(self):
                Joint.__init__(self)
                self.dofs = [0, 0, 0, 1, 1, 1]
                self.name = 'spherical-joint'
                



class RevoluteJoint(Joint):

        def __init__(self, axis):
                Joint.__init__(self)
                self.dofs[3 + axis] = 1
                self.name = 'revolute-joint'

class CylindricalJoint(Joint):

        def __init__(self, axis ):
                Joint.__init__(self)
                self.dofs[0 + axis] = 1
                self.dofs[3 + axis] = 1
                self.name = 'cylindrical-joint'

class PrismaticJoint(Joint):

        def __init__(self, axis):
                Joint.__init__(self)
                self.dofs[0 + axis] = 1
                self.name = 'prismatic-joint'
