//BPRD Assignment 11 by Casper Brock (casb@itu.dk)

sealed abstract class Expr
case class CstI(value: Int) extends Expr
case class Var(name: String) extends Expr              //Exercise 11.1 (ii)
case class Prim(op: String, e1: Expr, e2: Expr) extends Expr


object exercise11 {
def main(args: Array[String]) = { 
  println(eval0(Prim("*", CstI(7), Prim("-", CstI(10), Prim("+", CstI(5), CstI(6))))))  //7 * (10 - 6 + 5)
  println(eval0(Prim("*", CstI(7), Prim("<", CstI(10), Prim("+", CstI(5), CstI(6))))))  //7 * (10 < 6 + 5)
  println(eval1(Prim("+", Var("x"), Var("y")), Map("x" -> 1, "y" -> 2)))   //Testing 11.1 (ii)
}

def eval0(expr: Expr): Int = 
  expr match {
    case CstI(i) => i
    case Prim(op, e1, e2) =>
      val v1 = eval0(e1)
      val v2 = eval0(e2)
        op match {
          case "+" => v1 + v2
          case "-" => v1 - v2
          case "*" => v1 * v2
          case "/" if(v2 != 0) => v1 / v2               //Exercise 11.1 (i)
          case "<" => if(v1 < v2) 1 else 0              //Exercise 11.1 (i)
          case ">" => if(v1 > v2) 1 else 0              //Exercise 11.1 (i)
          case "&" => if(v1 != 0 && v2 != 0) 1 else 0   //Exercise 11.1 (i)
      }
}

def eval1(expr: Expr, env: Map[String, Int]): Int = 
  expr match {
    case CstI(i) => i
    case Var(x) => env get x match {
                           case Some(v) => v
                          }                 
    case Prim(op, e1, e2) =>
      val v1 = eval1(e1, env)
      val v2 = eval1(e2, env)
        op match {
          case "+" => v1 + v2
          case "-" => v1 - v2
          case "*" => v1 * v2
          case "/" if(v2 != 0) => v1 / v2               //Exercise 11.1
          case "<" => if(v1 < v2) 1 else 0              //Exercise 11.1
          case ">" => if(v1 > v2) 1 else 0              //Exercise 11.1
          case "&" => if(v1 != 0 && v2 != 0) 1 else 0   //Exercise 11.1
      }
}

def simplify(expr: Expr): Expr =                        //Exercise 11.2
  expr match {
    case Prim("+", CstI(0), e) => e
    case Prim("+", e, CstI(0)) => e
    case Prim("-", e, CstI(0)) => e
    case Prim("*", e, CstI(1)) => e
    case Prim("*", CstI(1), e) => e
    case Prim("*", e, CstI(0)) => CstI(0)
    case Prim("*", CstI(0), e) => CstI(0)
    case Prim("-", e1, e2) if(e1 == e2) => CstI(0)
    case Prim(op, e1, e2) => Prim(op, simplify(e1), simplify(e2)) 
    case _ => expr
}

def eval2(expr: Expr, env: Map[String, Int]): Option[Int] =  //Exercise 11.3
  expr match {
    case CstI(i) => Some(i)
    case Var(x) => if(env contains x)
                      env get x
                   else None                
    case Prim(op, e1, e2) =>
      val v1 = eval2(e1, env)
      val v2 = eval2(e2, env)
        v1 match {
           case None => None
           case Some(i1) => {
              v2 match {
                 case None => None
                 case Some(i2) => {
                      op match {
                         case "+" => Some(i1 + i2)
                         case "-" => Some(i1 - i2)
                         case "*" => Some(i1 * i2)
                         case "/" if(i2 != 0) => Some(i1 / i2)
                         case "/" if(i2 == 0) => None
                         case "<" => if(i1 < i2) Some(1) else Some(0)
                         case ">" => if(i1 > i2) Some(1) else Some(0)
                         case "&" => if(i1 != 0 && i2 != 0) Some(1) else Some(0)
      }
     }
    }
   }
  }
 }

def eval3(e: Expr, env: Map[String, Int]): Option[Int] = e match {
  case CstI(n) => Some(n)
  case Var(x) => env get x
  case Prim(op, e1, e2) =>
     for {
       v1 <- eval3(e1, env)
       v2 <- eval3(e2, env)
       res <- op match {
         case "+" => Some(v1 + v2)
         case "-" => Some(v1 - v2)
         case "*" => Some(v1 * v2)
         case "/" => if (v2 == 0) None else Some(v1 / v2)
         case ">" => Some(if (v1 > v2) 1 else 0)
         case "<" => Some(if (v1 < v2) 1 else 0)
         case "&" => Some(v1 & v2)
       }
     } yield res
   }
}