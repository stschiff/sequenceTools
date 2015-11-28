module OrderedZip (orderedZip) where

import Control.Concurrent (myThreadId)
import qualified Control.Concurrent.STM as STM
import Control.Concurrent.Async (withAsync, wait)
import Foreign.StablePtr (newStablePtr)
import Turtle

data Slot a = Empty | Value a | Done

orderedZip :: (a -> b -> Ordering) -> Shell a -> Shell b -> Shell (Maybe a, Maybe b)
orderedZip cmp sA sB = Shell _foldIOAB
  where
    _foldIOAB (FoldM stepAB beginAB doneAB) = do
        -- myThreadId >>= newStablePtr
        x0 <- beginAB

        tvar <- STM.atomically (STM.newTVar (Empty, Empty))

        let begin = return ()

        let stepA () a = STM.atomically (do
                (x, y) <- STM.readTVar tvar
                case x of
                    Empty -> STM.writeTVar tvar (Value a, y)
                    _     -> STM.retry )
        let doneA () = STM.atomically (do
                (x, y) <- STM.readTVar tvar
                case x of
                    Empty -> STM.writeTVar tvar (Done, y)
                    _     -> STM.retry )
        let foldA = FoldM stepA begin doneA

        let stepB () b = STM.atomically (do
                (x, y) <- STM.readTVar tvar
                case y of
                    Empty -> STM.writeTVar tvar (x, Value b)
                    _      -> STM.retry )
        let doneB () = STM.atomically (do
                (x, y) <- STM.readTVar tvar
                case y of
                    Empty -> STM.writeTVar tvar (x, Done)
                    _      -> STM.retry )
        let foldB = FoldM stepB begin doneB

        withAsync (foldIO sA foldA) (\asyncA -> do
            withAsync (foldIO sB foldB) (\asyncB -> do
                let loop x = do
                        y <- STM.atomically (do
                            z <- STM.readTVar tvar
                            case z of
                                (Value a, Value b) -> do
                                    case a `cmp` b of
                                        LT -> do
                                            STM.writeTVar tvar (Empty, Value b)
                                            return $ Just (Just a, Nothing)
                                        EQ -> do
                                            STM.writeTVar tvar (Empty, Empty)
                                            return $ Just (Just a, Just b)
                                        GT -> do
                                            STM.writeTVar tvar (Value a, Empty)
                                            return $ Just (Nothing, Just b)
                                (Done, Value b) -> do
                                    STM.writeTVar tvar (Done, Empty)
                                    return $ Just (Nothing, Just b)
                                (Value a, Done) -> do
                                    STM.writeTVar tvar (Empty, Done)
                                    return $ Just (Just a, Nothing)
                                (Done, Done) -> return Nothing
                                _         -> STM.retry )
                        case y of
                            Nothing -> return x
                            Just ab -> do
                                x' <- stepAB x ab
                                loop $! x'
                x' <- loop $! x0
                wait asyncA
                wait asyncB
                doneAB x' ) )

testOrderedZip :: IO ()
testOrderedZip = do
    let shellA = select [1, 3, 6, 7, 9, 10, 15]
    let shellB = select [(2, "hello"), (3, "you"), (7, "test"), (11, "hi hi"), (15, "asdsa"), (16, "aaa")]
    let cmp a b = a `compare` (fst b)
    view (orderedZip cmp shellA shellB)